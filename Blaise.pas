(*

Pascal mathematica unit by Srbsialv D. Nešić, 1986--2020,
srbislav.nesic@gmail.com

*)

unit Blaise;


{$DEFINE WIRTH}


interface

const
  mathfpu = {$IFDEF WIRTH} false {$ELSE} true {$ENDIF};

type
  numeric = extended;
  ordinal = int64;
  natural = uint64;

type
  numarray = array of numeric;
  ordarray = array of ordinal;

function frac(x: numeric): numeric; overload;
function trunc(x: numeric): ordinal; overload;

function pi: numeric;

function isnan(const x: numeric): boolean; overload;
function isinf(const x: numeric): boolean; overload;
function isfin(const x: numeric): boolean; overload;
function ispos0(const x: numeric): boolean; overload;
function isneg0(const x: numeric): boolean; overload;

procedure order(var a, b: ordinal); overload;
procedure order(var a, b, c: ordinal); overload;
function min(m, n: ordinal): ordinal; overload;
function max(m, n: ordinal): ordinal; overload;
function min(a, b, c: ordinal): ordinal; overload;
function mid(a, b, c: ordinal): ordinal; overload;
function max(a, b, c: ordinal): ordinal; overload;
function censor(l, n, r: ordinal): ordinal; overload;
function abs(n: ordinal): ordinal; overload;
function sign(n: ordinal; k: ordinal = 1): ordinal; overload;
function gcd(m, n: ordinal): ordinal; overload;
function lcm(m, n: ordinal): ordinal; overload;
function bitcount(n: ordinal): integer; overload;
function basediv(n: natural): natural; overload;
function multinv(o: natural): natural; overload;
function hiadd(var carry: natural; a, b: natural): natural; overload;
function himul(var carry: natural; a, b: natural): natural; overload;
function fractional(x: numeric; var n, d: ordinal;
  depth: integer = 8; eps: numeric = 1e-12): numeric; overload;

function fibnumber(n: ordinal): numeric; overload;

procedure order(var a, b: numeric); overload;
procedure order(var a, b, c: numeric); overload;
function min(x, y: numeric): numeric; overload;
function max(x, y: numeric): numeric; overload;
function min(a, b, c: numeric): numeric; overload;
function mid(a, b, c: numeric): numeric; overload;
function max(a, b, c: numeric): numeric; overload;
function istriangle(a, b, c: numeric): boolean; overload;
function isright(a, b, c: numeric): boolean; overload;
function censor(l, x, r: numeric): numeric; overload;
function abs(x: numeric): numeric; overload;
function sign(x: numeric; y: numeric = 1): numeric; overload;
function psign(x: numeric; y: numeric = 1): numeric; overload;
function nsign(x: numeric; y: numeric = 1): numeric; overload;
function zsign(x: numeric; y: numeric = 1): numeric; overload;
function recip(x: numeric; y: numeric = 1): numeric; overload;
function int(x: numeric): numeric; overload;
function nearint(x: numeric): numeric; overload;
function nearrec(x: numeric): numeric; overload;
function floor(x: numeric): ordinal; overload;
function ceil(x: numeric): ordinal; overload;
function round(x: numeric): ordinal; overload;
function modulo(x: numeric; y: numeric = 1): numeric; overload;
function wrap(x, l, r: numeric): numeric; overload;
function wrap(x, y: numeric): numeric; overload;
function wrap(x: numeric): numeric; overload;
function deg(rad: numeric = 1): numeric; overload;
function rad(deg: numeric = 1): numeric; overload;

function ilog2(var x: numeric): ordinal; // internal use

function sqr(x: numeric): numeric; overload;
function cub(x: numeric): numeric; overload;
function sqrt(x: numeric): numeric; overload;
function cbrt(x: numeric): numeric; overload;
function hypot(x, y: numeric): numeric; overload;
function compmod(x: numeric): numeric; overload;

function power(x: numeric; n: ordinal): numeric; overload;
function power(x, y: numeric): numeric; overload;

function log2(x: numeric): numeric; overload;
function ln(x: numeric): numeric; overload;
function log10(x: numeric): numeric; overload;
function log(x, y: numeric): numeric; overload;

function exp2(x: numeric): numeric; overload;
function exp(x: numeric): numeric; overload;
function exp10(x: numeric): numeric; overload;

function naplog(x: numeric; n: numeric = 10000000): numeric; overload;
function napexp(x: numeric; n: numeric = 10000000): numeric; overload;
function mflog(x: numeric): numeric; overload;
function mfexp(x: numeric): numeric; overload;
function lnone(x: numeric): numeric; overload;
function logit(x: numeric; a: numeric = 1): numeric; overload;
function logistic(x: numeric; a: numeric = 1): numeric; overload;

procedure sincos(x: numeric; var sin, cos: numeric); overload;
function sin(x: numeric): numeric; overload;
function cos(x: numeric): numeric; overload;
function tan(x: numeric): numeric; overload;
function cot(x: numeric): numeric; overload;
function sec(x: numeric): numeric; overload;
function csc(x: numeric): numeric; overload;
function hav(x: numeric): numeric; overload;
function crd(x: numeric): numeric; overload;
function sinc(x: numeric): numeric; overload;

function arg(x, y: numeric): numeric; overload;
function atan2(y, x: numeric): numeric; overload;
function arcsin(x: numeric): numeric; overload;
function arccos(x: numeric): numeric; overload;
function arctan(x: numeric): numeric; overload;
function arccot(x: numeric): numeric; overload;
function arccsc(x: numeric): numeric; overload;
function arcsec(x: numeric): numeric; overload;
function archav(x: numeric): numeric; overload;
function arccrd(x: numeric): numeric; overload;
function arctanloc(x, location: numeric): numeric; overload;

procedure sindcosd(x: numeric; var sin, cos: numeric); overload;
function sind(x: numeric): numeric; overload;
function cosd(x: numeric): numeric; overload;
function tand(x: numeric): numeric; overload;
function cotd(x: numeric): numeric; overload;
function argd(x, y: numeric): numeric; overload;
function arcsind(x: numeric): numeric; overload;
function arccosd(x: numeric): numeric; overload;
function arctand(x: numeric): numeric; overload;
function arccotd(x: numeric): numeric; overload;

procedure sinhcosh(x: numeric; var sinh, cosh: numeric); overload;
function sinh(x: numeric): numeric; overload;
function cosh(x: numeric): numeric; overload;
function tanh(x: numeric): numeric; overload;
function coth(x: numeric): numeric; overload;

function arcsinh(x: numeric): numeric; overload;
function arccosh(x: numeric): numeric; overload;
function arctanh(x: numeric): numeric; overload;
function arccoth(x: numeric): numeric; overload;

function hadd(x, y: numeric): numeric; overload;
function hsub(x, y: numeric): numeric; overload;
function hmean(x, y: numeric): numeric; overload;
function amean(x, y: numeric; t: numeric = 1/2): numeric; overload;
function gmean(x, y: numeric): numeric; overload;
function agm(x, y: numeric): numeric; overload;

function EllipticK(x: numeric): numeric; overload;

function PendulumAmplitudePeriod(length, inclination: numeric;
  gravity: numeric = 9.80665): numeric; overload;
function TimeDilation(v: numeric): numeric; overload;
function Gudermannian(x: numeric): numeric; overload;
function AntiCollision(n, k:integer): numeric; overload;
function Gaussian(x: numeric;
  a:numeric = 1; b:numeric = 0; c: numeric = 1; d: numeric = 0; p: numeric = 2): numeric; overload;

function mvhgdprob(const play, hits: array of integer): numeric; overload;
function hgdprob(balls, draws, play, hits: integer): numeric; overload;
function KenoProb(play, hits: integer): numeric; overload;

function midifreq(note: numeric): numeric;
function midinote(freq: numeric): numeric;
function pianofreq(key: numeric): numeric;
function pianokey(freq: numeric): numeric;

function Kahan(var s, o, c: numeric; a: numeric): boolean; overload;

type
  complex = packed record
    x, y: numeric;
  private
    function getabs: numeric;
    function getarg: numeric;
    function getreal: boolean;
    function getimag: boolean;
    procedure setreal(f: boolean = true);
    procedure setimag(f: boolean = true);
    function getintX: ordinal;
    function getintY: ordinal;
    property isreal: boolean read getreal write setreal;
    property isimag: boolean read getimag write setimag;
  public
    class operator implicit(x: numeric): complex; overload;
    class operator equal(u, v: complex): boolean;
    class operator notequal(u, v: complex): boolean;
    class operator positive(z: complex): complex;
    class operator negative(z: complex): complex;
    class operator logicalnot(z: complex): complex;
    class operator add(x: numeric; z: complex): complex; overload;
    class operator add(z: complex; x: numeric): complex; overload;
    class operator add(u, v: complex): complex; overload;
    class operator subtract(x: numeric; z: complex): complex; overload;
    class operator subtract(z: complex; x: numeric): complex; overload;
    class operator subtract(u, v: complex): complex; overload;
    class operator multiply(z: complex; x: numeric): complex; overload;
    class operator multiply(x: numeric; z: complex): complex; overload;
    class operator multiply(u, v: complex): complex; overload;
    class operator divide(z: complex; x: numeric): complex; overload;
    class operator divide(x: numeric; z: complex): complex; overload;
    class operator divide(u, v: complex): complex; overload;
    property ρ: numeric read getabs;
    property θ: numeric read getarg;
    property ϑ: numeric read getarg;
    property Χ: ordinal read getintX;
    property Υ: ordinal read getintY;
  end;

const
  ί: complex = (x: 0; y: 1);

function isnan(const z: complex): boolean; overload;
function isinf(const z: complex): boolean; overload;
function isfin(const z: complex): boolean; overload;
function xiy(x: numeric; y: numeric = 0): complex; overload;
function cis(θ: numeric): complex; overload;
function cis(θ: numeric; var x, y: numeric): complex; overload;
function taxicab(z: complex): numeric; overload;
function sqrabs(z: complex): numeric; overload;
function abs(z: complex): numeric; overload;
function arg(z: complex): numeric; overload;
function cisd(θ: numeric): complex; overload;
function argd(z: complex): numeric; overload;
function eis(z: complex): complex; overload;
function sign(z: complex): complex; overload;
function sign(x: numeric; z: complex): complex; overload;
function sqr(z: complex): complex; overload;
function sqrt(z: complex): complex; overload;
function compmod(z: complex): complex; overload;
function floor(z: complex): complex; overload;
function ceil(z: complex): complex; overload;
function round(z: complex): complex; overload;
function modulo(u, v: complex): complex; overload;

function fibnumber(x: numeric): complex; overload;

function dotcross(u, v: complex): complex; overload;
function dotprod(u, v: complex): numeric; overload;
function crossprod(u, v: complex): numeric; overload;
function intersection(u1, u2, v1, v2: complex): complex; overload;

function amean(u, v: complex; t: numeric = 1/2): complex; overload;
function gmean(u, v: complex): complex; overload;
function agm(u, v: complex): complex; overload;
function curve(u, p, v: complex; t: numeric = 1/2): complex; overload;
function bezier(u, p, q, v: complex; t: numeric = 1/2): complex;

function exp(z: complex): complex; overload;
function ln(z: complex): complex; overload;
function power(z: complex; n: ordinal): complex; overload;
function power(z: complex; x: numeric): complex; overload;
function power(u, v: complex): complex; overload;
function power(x: numeric; z: complex): complex; overload;
function log2(z: complex): complex; overload;
function exp2(z: complex): complex; overload;
function root(z: complex; n: ordinal): complex; overload;

procedure sinhcosh(z: complex; var sinh, cosh: complex); overload;
function sinh(z: complex): complex; overload;
function cosh(z: complex): complex; overload;
function tanh(z: complex): complex; overload;
function coth(z: complex): complex; overload;

procedure sincos(z: complex; var sin, cos: complex); overload;
function sin(z: complex): complex; overload;
function cos(z: complex): complex; overload;
function tan(z: complex): complex; overload;
function cot(z: complex): complex; overload;
function sec(z: complex): complex; overload;
function csc(z: complex): complex; overload;
function hav(z: complex): complex; overload;
function crd(z: complex): complex; overload;
function sinc(z: complex): complex; overload;

function arcsinh(z: complex): complex; overload;
function arccosh(z: complex): complex; overload;
function arctanh(z: complex): complex; overload;
function arccoth(z: complex): complex; overload;

function arcsin(z: complex): complex; overload;
function arccos(z: complex): complex; overload;
function arctan(z: complex): complex; overload;
function arccot(z: complex): complex; overload;
function arccsc(z: complex): complex; overload;
function arcsec(z: complex): complex; overload;
function archav(z: complex): complex; overload;
function arccrd(z: complex): complex; overload;
function arctanloc(z, location: complex): complex; overload;

function PendulumAmplitudePeriod(body: complex; gravity: numeric = 9.80665): numeric; overload;

function LambertW(x: numeric): numeric; overload;
function LambertW(z: complex): complex; overload;

function poly(x: numeric; const a: array of numeric): numeric; overload;
function poly(z: complex; const a: array of numeric): complex; overload;
function poly(z: complex; const a: array of complex): complex; overload;

function TropicalYear(year: numeric): numeric; overload;

procedure exchange(var a, b: ordinal); overload;
procedure exchange(var a, b: numeric); overload;
procedure exchange(var a, b: complex); overload;

function Gamma(x: numeric): numeric; overload;
function Gamma(z: complex): complex; overload;
function fact(x: numeric): numeric; overload;
function fact(z: complex): complex; overload;
function fact2(x: numeric): numeric; overload;
function fact2(z: complex): complex; overload;
function binom(n, k: ordinal): numeric; overload;
function binom(n, k: numeric): numeric; overload;
function binom(n, k: complex): complex; overload;
function multinom(const x: array of numeric): numeric; overload;
function multinom(const z: array of complex): complex; overload;
function Beta(u, v: ordinal): numeric; overload;
function Beta(u, v: numeric): numeric; overload;
function Beta(u, v: complex): complex; overload;

function fallfact(x, n: numeric): numeric; overload;
function risefact(x, n: numeric): numeric; overload;
function fallfact(z, n: complex): complex; overload;
function risefact(z, n: complex): complex; overload;

function BallVolume(dimension: numeric): numeric; overload;
function SphereArea(dimension: numeric): numeric; overload;
function Lissajous(t, a, b, u, v, d: numeric): complex; overload;
function Lissajous(t, a, b, u, v: numeric): complex; overload;
function Lissajous(t, u, v: numeric): complex; overload;

function hms(h: numeric; m: numeric = 0; s: numeric = 0): numeric; overload;
function dms(d: numeric; m: numeric = 0; s: numeric = 0): numeric; overload;
function GeoPos(lat, lon: numeric): complex; overload;
function GeoDist(a, b: complex; r: numeric = 6371008.8): numeric; overload;

function mandelbrot(c: complex; x: numeric = 2; m: integer = 1000): integer; overload;
function mandelbrot(c, p: complex; m: integer = 1000): integer; overload;

function superellipse(radius: complex; angle, shape, symmetry, u, v: numeric): complex; overload;
function superellipse(radius: complex; angle: numeric;
  shape: numeric = 2; symmetry: numeric = 4): complex; overload;
function superellipse(radius: numeric; angle: numeric;
  shape: numeric = 2; symmetry: numeric = 4): complex; overload;

type
  dot = packed record
    x, y: ordinal;
    class operator implicit(d: dot): complex; overload;
    class operator implicit(z: complex): dot; overload;
  end;

type
  rational = packed record
    p, q: ordinal;
  private
    function getx: numeric;
    procedure setx(r: numeric);
  public
    class operator implicit(x: numeric): rational; overload;
    class operator implicit(n: ordinal): rational; overload;
    class operator implicit(f: rational): numeric; overload;
    class operator implicit(f: rational): complex; overload;
    class operator equal(f, g: rational): boolean;
    class operator notequal(f, g: rational): boolean;
    class operator lessthan(f, g: rational): boolean;
    class operator greaterthan(f, g: rational): boolean;
    class operator lessthanorequal(f, g: rational): boolean;
    class operator greaterthanorequal(f, g: rational): boolean;
    class operator positive(f: rational): rational;
    class operator negative(f: rational): rational;
    class operator logicalnot(f: rational): rational;
    class operator multiply(f, g: rational): rational;
    class operator divide(f, g: rational): rational;
    class operator intdivide(f, g: rational): ordinal;
    class operator modulus(f, g: rational): rational;
    class operator add(f, g: rational): rational;
    class operator subtract(f, g: rational): rational;
    class operator inc(f: rational): rational; overload;
    class operator dec(f: rational): rational; overload;
    class operator trunc(f: rational): ordinal; overload;
    property x: numeric read getx write setx;
  end;

type ratarray = array of rational;

function isnan(const f: rational): boolean; overload;
function isinf(const f: rational): boolean; overload;
function isfin(const f: rational): boolean; overload;
function frc(p: ordinal; q: ordinal = 1): rational; overload;
function trunc(f: rational): ordinal; overload;
function frac(f: rational): rational; overload;
function floor(f: rational): ordinal; overload;
function ceil(f: rational): ordinal; overload;
function round(f: rational): ordinal; overload;
function sign(f: rational): ordinal; overload;
function abs(f: rational): rational; overload;
procedure order(var a, b: rational); overload;
procedure order(var a, b, c: rational); overload;
function min(f, g: rational): rational; overload;
function max(f, g: rational): rational; overload;
function min(a, b, c: rational): rational; overload;
function mid(a, b, c: rational): rational; overload;
function max(a, b, c: rational): rational; overload;
function censor(l, n, r: rational): rational; overload;
function gcd(f, g: rational): rational; overload;
function lcm(f, g: rational): rational; overload;
function sqr(f: rational): rational; overload;
function amean(f, g: rational): rational; overload;
function hadd(f, g: rational): rational; overload;
function hsub(f, g: rational): rational; overload;
function hmean(f, g: rational): rational; overload;
function contfrac(const v: array of ordinal; n: integer = -1): rational; overload;
procedure exchange(var a, b: rational); overload;
function fractional(x: numeric; var f: rational;
  depth: integer = 8; eps: numeric = 1e-12): numeric; overload;

function simplify(f: rational; depth: integer = 8): rational; overload;
function crossprod(f, g: rational): ordinal; overload;
procedure fareyseq(n: integer; var f: ratarray);
function arctan(f: rational): numeric; overload;


type
  quaternion =  packed record
  private
    function getisscalar: boolean;
    procedure setisscalar(f: boolean);
    function getisvector: boolean;
    procedure setisvector(f: boolean);
    function geta: complex;
    procedure seta(z: complex);
    function getisa: boolean;
    function getb: complex;
    procedure setb(z: complex);
    function getisb: boolean;
    function getc: complex;
    procedure setc(z: complex);
    function getisc: boolean;
    property isscalar: boolean read getisscalar write setisscalar;
    property isvector: boolean read getisvector write setisvector;
    property isa: boolean read getisa;
    property isb: boolean read getisb;
    property isc: boolean read getisc;
  public
    property a: complex read geta write seta;
    property b: complex read getb write setb;
    property c: complex read getc write setc;
    class operator implicit(x: numeric): quaternion; overload;
    class operator implicit(z: complex): quaternion; overload;
    class operator positive(q: quaternion): quaternion; overload;
    class operator negative(q: quaternion): quaternion; overload;
    class operator logicalnot(q: quaternion): quaternion; overload;
    class operator equal(p, q: quaternion): boolean; overload;
    class operator notequal(p, q: quaternion): boolean; overload;
    class operator add(p, q: quaternion): quaternion; overload;
    class operator add(q: quaternion; x: numeric): quaternion; overload;
    class operator add(x: numeric; q: quaternion): quaternion; overload;
    class operator subtract(p, q: quaternion): quaternion; overload;
    class operator subtract(q: quaternion; x: numeric): quaternion; overload;
    class operator subtract(x: numeric; q: quaternion): quaternion; overload;
    class operator multiply(x: numeric; q: quaternion): quaternion; overload;
    class operator multiply(q: quaternion; x: numeric): quaternion; overload;
    class operator multiply(p, q: quaternion): quaternion; overload; // not commutative
    class operator divide(q: quaternion; x: numeric): quaternion; overload;
    class operator divide(x: numeric; q: quaternion): quaternion; overload;
    class operator divide(p, q: quaternion): quaternion; overload;
    class operator bitwiseand(p, q: quaternion): quaternion; overload; // not commutative
    class operator leftshift(p, q: quaternion): quaternion; overload;
    class operator rightshift(p, q: quaternion): quaternion; overload;
  case integer of
    1: (w, x, y, z: numeric); // w + x·i + y·j + z·k, i² = j² = k² = ijk = -1
    2: (s: numeric; v: packed array [1..3] of numeric);
    3: (o: packed array [0..3] of numeric);
    4: (f, g: complex); // f + g·j
  end;

const
  і: quaternion = (w: 0; x: 1; y: 0; z: 0);
  ϳ: quaternion = (w: 0; x: 0; y: 1; z: 0);
  ƙ: quaternion = (w: 0; x: 0; y: 0; z: 1);

function isnan(const q: quaternion): boolean; overload;
function isinf(const q: quaternion): boolean; overload;
function isfin(const q: quaternion): boolean; overload;
function qtn(w: numeric; x: numeric = 0; y: numeric = 0; z: numeric = 0): quaternion; overload;
function qtn(u, v: complex): quaternion; overload;
function xyz(x, y, z: numeric): quaternion; overload;
function sqrabs(q: quaternion): numeric; overload;
function abs(q: quaternion): numeric; overload;
function arg(q: quaternion): quaternion; overload;
function sign(q: quaternion): quaternion; overload;
function versor(q: quaternion): quaternion; overload; // v² = -1
function polar(q: quaternion): quaternion; overload;
function polar(x: numeric; y: numeric = 0; z: numeric = 0): quaternion; overload; // φ, θ, ψ
function rotation(q: quaternion): quaternion; overload;
function rotation(x: numeric; y: numeric = 0; z: numeric = 0): quaternion; overload;
function angles(q: quaternion): quaternion; overload;
function vectorpart(q: quaternion): quaternion; overload;
function scalarpart(q: quaternion): quaternion; overload;
function sqr(q: quaternion): quaternion; overload;
function sqrt(q: quaternion): quaternion; overload;
function exp(q: quaternion): quaternion; overload;
function ln(q: quaternion): quaternion; overload;

function power(q: quaternion; x: numeric): quaternion; overload;
function power(q: quaternion; z: complex): quaternion; overload;
function power(p, q: quaternion): quaternion; overload;
function amean(p, q: quaternion; t: numeric = 1/2): quaternion; overload;
function dotcross(p, q: quaternion): quaternion; overload;
function dotprod(p, q: quaternion): numeric; overload;
function crossprod(p, q: quaternion): quaternion; overload;
function intersection(p, q: quaternion): complex; overload;
procedure exchange(var a, b: quaternion); overload;

function sin(q: quaternion): quaternion; overload;
function cos(q: quaternion): quaternion; overload;
function tan(q: quaternion): quaternion; overload;
function cot(q: quaternion): quaternion; overload;
function arcsin(q: quaternion): quaternion; overload;
function arccos(q: quaternion): quaternion; overload;
function arctan(q: quaternion): quaternion; overload;
function arccot(q: quaternion): quaternion; overload;

function sinh(q: quaternion): quaternion; overload;
function cosh(q: quaternion): quaternion; overload;
function tanh(q: quaternion): quaternion; overload;
function coth(q: quaternion): quaternion; overload;
function arcsinh(q: quaternion): quaternion; overload;
function arccosh(q: quaternion): quaternion; overload;
function arctanh(q: quaternion): quaternion; overload;
function arccoth(q: quaternion): quaternion; overload;

function Gamma(q: quaternion): quaternion; overload;

function mandelbrot(c: quaternion; x: numeric = 2; m: integer = 1000): integer; overload;
function yprtoqtn(q: quaternion): quaternion; overload;
function yprtoqtn(yaw, pitch, roll: numeric): quaternion; overload;
function qtntoypr(q: quaternion): quaternion; overload;
function qtntoypr(q: quaternion; var yaw, pitch, roll: numeric): quaternion; overload;

procedure entropy(var seed: natural); overload;
function mmixrng(var seed: natural): natural; overload;
function ordrandom(var seed: natural; n: ordinal): ordinal; overload;
function ordrange(var seed: natural; n1, n2: ordinal): ordinal; overload;
function unidev(var seed: natural; x: numeric = 1): numeric; overload;
function unidev(var seed: natural; a, b: numeric): numeric; overload;
function expdev(var seed: natural; lambda: numeric = 1): numeric; overload;
function geomdev(var seed: natural; p: numeric): ordinal; overload;
function polardev(var seed: natural; r: numeric): complex; overload;
function gaussdev(var seed: natural): complex; overload;
function normaldev(var seed: natural; mu, sigma: complex; angle: numeric = 0): complex; overload;
function erlangdev(var seed: natural; k: integer; lambda: numeric = 1): numeric; overload;
function chi2dev(var seed: natural; nu: numeric): numeric; overload;
function studenttdev(var seed: natural; nu: numeric): numeric; overload;
function gammadev(var seed: natural; a: numeric): numeric; overload;
function poissondev(var seed: natural; lambda: numeric): ordinal; overload;
function toss(var seed: natural; p: numeric = 1/2): boolean; overload;
function berndev(var seed: natural; n, k: ordinal): boolean; overload;
function bindev(var seed: natural; n: ordinal; p: numeric = 1/2): ordinal; overload;
function benforddev(var seed: natural; m, n: ordinal): ordinal; overload;
procedure ordstir(var seed: natural; var x: array of ordinal; n: integer = -1); overload;
procedure ordfill(var seed: natural; var x: array of ordinal; ordered: boolean = false); overload;
procedure ordsample(var seed: natural; var x: array of ordinal; n: integer); overload;

function mmixrng: natural; overload;
function ordrandom(n: ordinal): ordinal; overload;
function ordrange(n1, n2: ordinal): ordinal; overload;
function unidev(x: numeric = 1): numeric; overload;
function unidev(a, b: numeric): numeric; overload;
function expdev(lambda: numeric = 1): numeric; overload;
function geomdev(p: numeric): ordinal; overload;
function polardev(r: numeric = 1): complex; overload;
function gaussdev: complex; overload;
function normaldev(mu, sigma: complex; angle: numeric = 0): complex; overload;
function erlangdev(k: integer; lambda: numeric = 1): numeric; overload;
function chi2dev(nu: numeric): numeric; overload;
function studenttdev(nu: numeric): numeric; overload;
function gammadev(a: numeric): numeric; overload;
function poissondev(lambda: numeric): ordinal; overload;
function toss(p: numeric = 1/2): boolean; overload;
function berndev(n, k: ordinal): boolean; overload;
function bindev(n: ordinal; p: numeric = 1/2): ordinal; overload;
function benforddev(m, n: ordinal): ordinal; overload;
procedure ordstir(var x: array of ordinal; n: integer = -1); overload;
procedure ordfill(var x: array of ordinal; ordered: boolean = false); overload;
procedure ordsample(var x: array of ordinal; n: integer); overload;

procedure ordsort(var x: array of ordinal; lpart, rpart: integer); overload;
procedure ordsort(var x: array of ordinal); overload;
procedure numsort(var x: array of numeric; lpart, rpart: integer); overload;
procedure numsort(var x: array of numeric); overload;



implementation


// For implementation send mail to author


initialization


init;


finalization


done;


end.
