param   t1 = "Small Monte Carlo by Scott Prahl (https://omlc.org)";
param   t2 = "1 W/cm^2 Uniform Illumination of Semi-Infinite Medium";

use Math;
use Random;
use ChplConfig;
use GPU;

param BINS = 101;

config const useGpu = if CHPL_LOCALE_MODEL=="gpu" then here.gpus.size>0
                                                  else false;

config const   numPhotons = 100000;

record photon {
  var x,y,z,u,v,w,weight: real;
}

require "curand_kernel.h";
require "kernel_rng.h";
extern type curandState_t;

pragma "codegen for CPU and GPU"
extern proc rng_init(idx, seed, ref state: curandState_t): void;
pragma "codegen for CPU and GPU"
extern proc rng_get(ref state: curandState_t): real(32);

param RAND_MAX = max(int(32));
proc rand(ref rng: curandState_t) {
  const ret = (rng_get(rng)*RAND_MAX):int(32);
  return ret;
}

proc ref photon.launch(rs) /* Start the photon */
{
  x = 0.0; y = 0.0; z = 0.0;
  u = 0.0; v = 0.0; w = 1.0;
  weight = 1.0 - rs;
}

proc ref photon.bounce (n, crit_angle, ref rd) /* Interact with top surface */
{
  var t, temp, temp1,rf: real;
  w = -w;
  z = -z;
  if w <= crit_angle then return;		/* total internal reflection */

  t       = sqrt(1.0-n*n*(1.0-w*w));            /* cos of exit angle */
  temp1   = (w - n*t)/(w + n*t);
  temp    = (t - n*w)/(t + n*w);
  rf      = (temp1*temp1+temp*temp)/2.0;	/* Fresnel reflection */
  rd     += (1.0-rf) * weight;
  weight -= (1.0-rf) * weight;
}

proc ref photon.move(n, crit_angle, ref rng, ref rd) /* move to next scattering or absorption event */
{
  var d = -log((rand(rng)+1.0)/(RAND_MAX+1.0)); // TODO
  /*var d = 1;*/
  x += d * u;
  y += d * v;
  z += d * w;
  if z<=0 then bounce(n, crit_angle, rd);
}

proc ref photon.absorb (ref heat, bins_per_mfp, albedo, ref rng, ref bit) /* Absorb light in the medium */
{
  var bin=(z*bins_per_mfp):int;

  if bin >= BINS then bin = BINS-1;
  gpuAtomicAdd(heat[bin], (1.0-albedo)*weight);
  weight *= albedo;
  if weight < 0.001 { /* Roulette */
    bit -= weight;
    if rand(rng) > 0.1*RAND_MAX then weight = 0; else weight /= 0.1; //TODO
    bit += weight;
  }
}

proc ref photon.scatter(g, ref rng) /* Scatter photon and establish new direction */
{
  var x1, x2, x3, t, mu: real;

  while true {						/*new direction*/
    x1=2.0*rand(rng)/RAND_MAX - 1.0;  // TODO
    x2=2.0*rand(rng)/RAND_MAX - 1.0;  // TODO
    x3=x1*x1+x2*x2;
    if x3<=1 then break;
  }
  if (g==0) {  /* isotropic */
    u = 2.0 * x3 -1.0;
    v = x1 * sqrt((1-u*u)/x3);
    w = x2 * sqrt((1-u*u)/x3);
    return;
  }

  mu = (1-g*g)/(1-g+2.0*g*rand(rng)/RAND_MAX);  // TODO
  mu = (1 + g*g-mu*mu)/2.0/g;
  if abs(w) < 0.9 {
    t = mu * u + sqrt((1-mu*mu)/(1-w*w)/x3) * (x1*u*w-x2*v);
    v = mu * v + sqrt((1-mu*mu)/(1-w*w)/x3) * (x1*v*w+x2*u);
    w = mu * w - sqrt((1-mu*mu)*(1-w*w)/x3) * x1;
  } else {
    t = mu * u + sqrt((1-mu*mu)/(1-v*v)/x3) * (x1*u*v + x2*w);
    w = mu * w + sqrt((1-mu*mu)/(1-v*v)/x3) * (x1*v*w - x2*u);
    v = mu * v - sqrt((1-mu*mu)*(1-v*v)/x3) * x1;
  }
  u = t;
}

proc main ()
{
  const mu_a = 5.0;			/* Absorption Coefficient in 1/cm */
  const mu_s = 95.0;			/* Scattering Coefficient in 1/cm */
  const g = 0.5;				/* Scattering Anisotropy -1<=g<=1 */
  const n = 1.5;				/* Index of refraction of medium */
  const microns_per_bin = 20;/* Thickness of one bin layer */
  const albedo = mu_s / (mu_s + mu_a);
  const rs = (n-1.0)*(n-1.0)/(n+1.0)/(n+1.0);	/* specular reflection */
  const crit_angle = sqrt(1.0-1.0/n/n);			/* cos of critical angle */
  const bins_per_mfp = 1e4/microns_per_bin/(mu_a+mu_s);

  writef("%s\n%s\n\nScattering = %8.3dr/cm\nAbsorption = %8.3dr/cm\n",t1,t2,mu_s,mu_a);
  writef("Anisotropy = %8.3dr\nRefr Index = %8.3dr\nPhotons    = %8i\n",g,n,numPhotons);

  const targetLoc = if useGpu then here.gpus[0] else here;

  on targetLoc {
    var heat: [0..<BINS] real;
    var rds, bits: [0..<numPhotons] real;

    @gpu.assertEligible
    forall (idx, rd, bit) in zip(0..<numPhotons, rds, bits) {
      var rng: curandState_t;
      rng_init(1, idx, rng);

      var p: photon;

      p.launch (rs);
      while p.weight > 0 {
        p.move (n, crit_angle, rng, rd);
        p.absorb (heat, bins_per_mfp, albedo, rng, bit);
        p.scatter (g, rng);
      }
    }

    const rd = + reduce rds;
    const bit = + reduce bits;

    writef("\n\nSpecular Refl      = %10.dr\nBackscattered Refl = %10.5dr",rs,rd/(bit+numPhotons));
    writef("\n\n Depth         Heat\n[microns]     [W/cm^3]\n");

    for i in 0..<BINS-1 {
      writef("%6.0dr    %12.5dr\n",i*microns_per_bin, heat[i]/microns_per_bin*1e4/(bit+numPhotons));
    }
    writef(" extra    %12.5dr\n",heat[BINS-1]/(bit+numPhotons));
  }
}
