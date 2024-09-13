param   t1 = "Small Monte Carlo by Scott Prahl (https://omlc.org)";
param   t2 = "1 W/cm^2 Uniform Illumination of Semi-Infinite Medium";

use Math;
use Random;
param BINS = 101;

const mu_a = 5.0;			/* Absorption Coefficient in 1/cm */
const mu_s = 95.0;			/* Scattering Coefficient in 1/cm */
const g = 0.5;				/* Scattering Anisotropy -1<=g<=1 */
const n = 1.5;				/* Index of refraction of medium */
const microns_per_bin = 20;/* Thickness of one bin layer */
var   i, photons = 100000;
var x,y,z,u,v,w,weight: real;
var rs, rd, bit, albedo, crit_angle, bins_per_mfp: real;
var heat: [0..<BINS] real;

/* to emulate C's rand() */
var random = new randomStream(int, seed=1);
extern const RAND_MAX: int(32);
proc rand() {
  return random.next(0, RAND_MAX);
}

proc launch() /* Start the photon */
{
  x = 0.0; y = 0.0; z = 0.0;
  u = 0.0; v = 0.0; w = 1.0;
  weight = 1.0 - rs;
}

proc bounce () /* Interact with top surface */
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

proc move() /* move to next scattering or absorption event */
{
  var d = -log((rand()+1.0)/(RAND_MAX+1.0));
  x += d * u;
  y += d * v;
  z += d * w;
  if z<=0 then bounce();
}

proc absorb () /* Absorb light in the medium */
{
  var bin=(z*bins_per_mfp):int;

  if bin >= BINS then bin = BINS-1;
  heat[bin] += (1.0-albedo)*weight;
  weight *= albedo;
  if weight < 0.001 { /* Roulette */
    bit -= weight;
    if rand() > 0.1*RAND_MAX then weight = 0; else weight /= 0.1;
    bit += weight;
  }
}

proc scatter() /* Scatter photon and establish new direction */
{
  var x1, x2, x3, t, mu: real;

  while true {						/*new direction*/
    x1=2.0*rand()/RAND_MAX - 1.0;
    x2=2.0*rand()/RAND_MAX - 1.0;
    x3=x1*x1+x2*x2;
    if x3<=1 then break;
  }
  if (g==0) {  /* isotropic */
    u = 2.0 * x3 -1.0;
    v = x1 * sqrt((1-u*u)/x3);
    w = x2 * sqrt((1-u*u)/x3);
    return;
  }

  mu = (1-g*g)/(1-g+2.0*g*rand()/RAND_MAX);
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

proc print_results() /* Print the results */
{
  writef("%s\n%s\n\nScattering = %8.3dr/cm\nAbsorption = %8.3dr/cm\n",t1,t2,mu_s,mu_a);
  writef("Anisotropy = %8.3dr\nRefr Index = %8.3dr\nPhotons    = %8i",g,n,photons);
  writef("\n\nSpecular Refl      = %10.dr\nBackscattered Refl = %10.5dr",rs,rd/(bit+photons));
  writef("\n\n Depth         Heat\n[microns]     [W/cm^3]\n");

  for i in 0..<BINS-1 {
    writef("%6.0dr    %12.5dr\n",i*microns_per_bin, heat[i]/microns_per_bin*1e4/(bit+photons));
  }
  writef(" extra    %12.5dr\n",heat[BINS-1]/(bit+photons));
}

proc main ()
{
  writeln(RAND_MAX);
  albedo = mu_s / (mu_s + mu_a);
  rs = (n-1.0)*(n-1.0)/(n+1.0)/(n+1.0);	/* specular reflection */
  crit_angle = sqrt(1.0-1.0/n/n);			/* cos of critical angle */
  bins_per_mfp = 1e4/microns_per_bin/(mu_a+mu_s);

  for i in 1..photons {
    launch ();
    while weight > 0 {
      move ();
      absorb ();
      scatter ();
    }
  }	
  print_results();
  return 0;
}
