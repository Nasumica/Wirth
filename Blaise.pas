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



type

  zxfloat = packed record

    a, e, d, c, b: byte;

  private

    function getx: numeric;
    procedure setx(r: numeric);
    function getm: numeric;
    function gete: integer;
  public

    class operator implicit(z: zxfloat): numeric; overload;

    class operator implicit(x: numeric): zxfloat; overload;

    class operator equal(f, g: zxfloat): boolean;

    class operator notequal(f, g: zxfloat): boolean;
    class operator lessthan(f, g: zxfloat): boolean;
    class operator greaterthan(f, g: zxfloat): boolean;
    class operator lessthanorequal(f, g: zxfloat): boolean;
    class operator greaterthanorequal(f, g: zxfloat): boolean;
    class operator positive(f: zxfloat): zxfloat;
    class operator negative(f: zxfloat): zxfloat;
    class operator add(f, g: zxfloat): zxfloat;
    class operator subtract(f, g: zxfloat): zxfloat;
    class operator multiply(f, g: zxfloat): zxfloat;
    class operator divide(f, g: zxfloat): zxfloat;
    property x: numeric read getx write setx;
    property mantissa: numeric read getm;
    property exponent: integer read gete;

  end;


function zxf(n: natural): zxfloat; overload;

function zxf(a, e, d, c, b: byte): zxfloat; overload;


const

  zxzero: zxfloat = (a:0; e: 0; d: 0; c: 0; b: 0);

  zxeps: zxfloat = (a:1; e: 0; d: 0; c: 0; b: 0); // 2^-128

  zxmin: zxfloat = (a:255; e: 255; d: 255; c: 255; b: 255);

  zxmax: zxfloat = (a:255; e: 127; d: 255; c: 255; b: 255);


type

  zxcomplex = packed record

    x, y: zxfloat;
    class operator implicit(x: numeric): zxcomplex; overload;

    class operator implicit(c: zxcomplex): complex; overload;

    class operator implicit(z: complex): zxcomplex; overload;

  end;


implementation


function frac(x: numeric): numeric;
begin
  result := system.frac(x);
end;

function trunc(x: numeric): ordinal;
begin
  result := system.trunc(x);
end;


var
  π: numeric = 0; // 3.1415926535897932384626433832795
  τ: numeric = 0; // 2π
  ℮: numeric = 0; // 2.7182818284590452353602874713527
  ln2: numeric = 0; // 0.69314718055994530941723212145818
  lg10: numeric = 0; // 3.3219280948873623478703194294894

var
  sqrt2, sqrth, sqrt3, sqrt5, sqrtπ, φ, nateps: numeric;

const
  ∞ = 1/0;
  Ø = 0/0;
const
  inf: numeric = ∞;
  nan: numeric = Ø;
const
  complexinf: complex = (x: ∞; y: ∞);
  complexnan: complex = (x: Ø; y: Ø);
const
  ratinf: rational = (p: 1; q: 0);
  ratnan: rational = (p: 0; q: 0);
const
  qtninf: quaternion = (w: ∞; x: ∞; y: ∞; z: ∞);
  qtnnan: quaternion = (w: Ø; x: Ø; y: Ø; z: Ø);


function bitequal(const x, y: array of byte): boolean; overload;
var
  i: integer;
begin
  result := length(x) = length(y);  i := high(x);
  while result and (i >= low(x)) do begin
    result := x[i] = y[i];
    i := i - 1;
  end;
end;

function bitequal(const x, y: numeric): boolean; overload;
type
  bytes = packed array [1 .. sizeof(numeric)] of byte;
var
  u: bytes absolute x;
  v: bytes absolute y;
begin
  result := bitequal(u, v);
end;

function isnan(const x: numeric): boolean; overload;
begin
  result := bitequal(x, nan) or bitequal(x, -nan);
end;

function isinf(const x: numeric): boolean; overload;
begin
  result := bitequal(x, inf) or bitequal(x, -inf);
end;

function isfin(const x: numeric): boolean; overload;
begin
  result := not (isnan(x) or isinf(x));
end;

function isnan(const z: complex): boolean; overload;
begin
  result := isnan(z.x) or isnan(z.y);
end;

function isinf(const z: complex): boolean; overload;
begin
  result := isinf(z.x) or isinf(z.y);
end;

function isfin(const z: complex): boolean; overload;
begin
  result := not (isnan(z) or isinf(z));
end;

function isnan(const q: quaternion): boolean; overload;
begin
  result := isnan(q.f) or isnan(q.g);
end;

function isinf(const q: quaternion): boolean; overload;
begin
  result := isinf(q.f) or isinf(q.g);
end;

function isfin(const q: quaternion): boolean; overload;
begin
  result := not (isnan(q) or isinf(q));
end;

function isnan(const f: rational): boolean; overload;
begin
  result := (f.q = 0) and (f.p = 0);
end;

function isinf(const f: rational): boolean; overload;
begin
  result := (f.q = 0) and (f.p <> 0);
end;

function isfin(const f: rational): boolean; overload;
begin
  result := f.q <> 0;
end;

function ispos0(const x: numeric): boolean; overload;
const y: numeric = 0 / 1;
begin
  result := bitequal(x, y);
end;

function isneg0(const x: numeric): boolean; overload;
const y: numeric = 0 / -1;
begin
  result := bitequal(x, y);
end;

function Kahan(var s, o, c: numeric; a: numeric): boolean;
begin
  o := s;
  a := a - c;
  s := s + a;
  c := s - o;
  c := c - a;
  result := s = o;
end;

{$IFDEF WIRTH}

function pi: numeric;
var
  n: integer;
  a, q, r, c: numeric;
begin
  if π = 0 then begin
    n := -8;  q := 16;  c := 0;
    repeat
      n := n + 8;  q := q / 16;
      a := 4/(n + 1) - 2/(n + 4) - 1/(n + 5) - 1/(n + 6);
    until Kahan(π, r, c, a * q);
  end;
  result := π;
end;

{$ELSE}

function pi: numeric;
asm
  FLDPI; FWAIT
end;

{$ENDIF}


var
  tinies: array of numeric; // 1/2^2^n, n = 0, 1, 2, ...
  expscale: ordinal; // max 2^2^n

procedure calctinies;
var
  s: numeric;
begin
  expscale := 1; s := 1/2;
  setlength(tinies, 0);
  while s > 0 do begin
    expscale := expscale * 2;
    setlength(tinies, length(tinies) + 1);
    tinies[high(tinies)] := s;
    s := sqr(s);
  end;
  expscale := expscale div 2;
end;

function ilog2(var x: numeric): ordinal; // internal use
var
  n: ordinal;
  k: integer;
  p: numeric;
begin
  result := 0;
  k := high(tinies);
  n := expscale;
  if x > 1 then begin
    while (x <> 1) and (k > low(tinies)) do begin
      n := n div 2;
      k := k - 1;
      p := 1/tinies[k];
      while x >= p do begin
        result := result + n;
        x := x / p;
      end;
    end;
    if x > sqrt2 then begin
      result := result + 1;
      x := x / 2;
    end;
  end else begin
  if x < 1 then
    while (x <> 1) and (k >= low(tinies)) do begin
      p := tinies[k];
      while x <= p do begin
        result := result - n;
        x := x / p;
      end;
      n := n div 2;
      k := k - 1;
    end;
    if x <= sqrth then begin
      result := result - 1;
      x := x * 2;
    end;
  end;
end;


function min(m, n: ordinal): ordinal;
begin
  if m < n then result := m else result := n;
end;

function max(m, n: ordinal): ordinal;
begin
  if m > n then result := m else result := n;
end;

procedure order(var a, b: ordinal); overload;
begin
  if a > b then exchange(a, b);
end;

procedure order(var a, b, c: ordinal); overload;
begin
  order(a, b);  order(a, c);  order(b, c);
end;

function min(a, b, c: ordinal): ordinal;
begin
  order(a, b, c);  result := a;
end;

function mid(a, b, c: ordinal): ordinal;
begin
  order(a, b, c);  result := b;
end;

function max(a, b, c: ordinal): ordinal;
begin
  order(a, b, c);  result := c;
end;

function censor(l, n, r: ordinal): ordinal;
begin
  if n < l then result := l else
  if n > r then result := r else
    result := n;
end;

function abs(n: ordinal): ordinal;
begin
  if n < 0 then result := -n else result := n;
end;

function sign(n: ordinal; k: ordinal = 1): ordinal;
begin
  if n < 0 then result := -k else
  if n > 0 then result :=  k else
    result := 0;
end;

function binom(n, k: ordinal): numeric;
var
  i: integer;
begin
  k := min(k, n - k);  result := 1;
  for i := 1 to k do result := result * (n - i + 1) / i;
end;

function gcd(m, n: ordinal): ordinal;
begin
  result := m;
  while n <> 0 do begin
    m := n;
    n := result mod n;
    result := m;
  end;
end;

function lcm(m, n: ordinal): ordinal;
begin
  if (m = 0) and (n = 0)
    then result := 0
    else result := m div gcd(m, n) * n;
end;

function bitcount(n: ordinal): integer;
begin
  result := 0;
  while n <> 0 do begin
    result := result + 1;
    n := n and (n - 1);
  end;
end;

function basediv(n: natural): natural; // 2⁶⁴ div n
begin
  if n > 1
    then result := (not n + 1) div n + 1
    else result := 0;
end;

function multinv(o: natural): natural; // o * result = 1  mod 2⁶⁴
var
  b, m: natural;
begin
  result := 0;
  if odd(o) then begin
    m := 0;  b := 1;
    while b <> 0 do begin
      m := m + b;
      if o * result and m <> 1 then
        result := result + b;
      b := m + 1;
    end;
  end;
end;

function hiadd(var carry: natural; a, b: natural): natural; // a + b = result + carry * 2⁶⁴
begin
  result := a + b;
  if (result < a) or (result < b)
    then carry := 1
    else carry := 0;
end;

function himul(var carry: natural; a, b: natural): natural; // a * b = result + carry * 2⁶⁴
const
  s = sizeof(natural) * 4;
var
  c, d: natural;
begin
  result := a * b;
  c := (a shl s) shr s;  a := a shr s;
  d := (b shl s) shr s;  b := b shr s;
  carry := (c * d) shr s;
  carry := hiadd(d, carry, a * d);
  carry := hiadd(c, carry, b * c);
  carry := carry shr s + (c + d) shl s + a * b;
end;

function fractional(x: numeric; var n, d: ordinal; depth: integer; eps: numeric): numeric;
const
  max = 16;
var
  f: array [0..max] of ordinal;
  g: ordinal;
  i, j: integer;
  a: numeric;
begin
  eps := abs(eps);  depth := min(depth, max);
  a := x;  n := round(a);  d := 1;
  i := 0;  f[i] := n;  result := n;
  while (abs(result - x) > eps) and (i < depth) do begin
    a := 1/(a - f[i]);  i := i + 1;  f[i] := round(a);
    n := 1;  d := 0;
    for j := i downto 0 do begin
      g := n;  n := f[j] * n + d;  d := g;
    end;
    result := n/d;
  end;
  if d < 0 then begin
    d := -d;  n := -n;
  end;
end;

function min(x, y: numeric): numeric;
begin
  if x < y then result := x else result := y;
end;

function max(x, y: numeric): numeric;
begin
  if x > y then result := x else result := y;
end;

procedure order(var a, b: numeric); overload;
begin
  if a > b then exchange(a, b);
end;

procedure order(var a, b, c: numeric); overload;
begin
  order(a, b);  order(a, c);  order(b, c);
end;

function istriangle(a, b, c: numeric): boolean;
begin
  a := abs(a);  b := abs(b);  c := abs(c);  order(a, b, c);
  result := ((a + b) > c) and (a > 0);
end;

function isright(a, b, c: numeric): boolean;
begin
  a := sqr(a);  b := sqr(b);  c := sqr(c);  order(a, b, c);
  result := ((a + b) = c) and (a > 0);
end;

function min(a, b, c: numeric): numeric;
begin
  order(a, b, c);  result := a;
end;

function mid(a, b, c: numeric): numeric;
begin
  order(a, b, c);  result := b;
end;

function max(a, b, c: numeric): numeric;
begin
  order(a, b, c);  result := c;
end;

function censor(l, x, r: numeric): numeric;
begin
  if x < l then result := l else
  if x > r then result := r else
    result := x;
end;

function abs(x: numeric): numeric;
begin
  if x < 0 then result := -x else result := x;
end;

function sign(x: numeric; y: numeric = 1): numeric;
begin
  if x < 0 then result := -y else
  if x > 0 then result :=  y else
    result := 0;
end;

function psign(x: numeric; y: numeric = 1): numeric;
begin
  if x < 0 then result := -y else result := y;
end;

function nsign(x: numeric; y: numeric = 1): numeric;
begin
  result := psign(-x, -y);
end;

function zsign(x: numeric; y: numeric = 1): numeric;
begin
  if isnan(x) or isnan(y) then result :=  0 else
  if (x < 0) or isneg0(x) then result := -y else
  if (x > 0) or ispos0(x) then result :=  y else
    result := 0;
end;

function recip(x: numeric; y: numeric = 1): numeric; overload;
begin
  result := zsign(x) * zsign(y);
  if result = 0 then result := Ø else begin
    x := abs(x);  y := abs(y);
    if x = 0 then begin
      if (y = 0) or (y = ∞)
        then result := Ø
        else result := sign(result, ∞);
    end else
    if x = ∞ then begin
      if (y = 0) or (y = ∞)
        then result := Ø
        else result := sign(result, 0)
    end else result := sign(result, y/x);
  end;
end;

function int(x: numeric): numeric;
begin
  result := x - frac(x);
end;

function nearint(x: numeric): numeric;
begin
  if isfin(x) then begin
    x := x + 1/2;  result := int(x);
    if result > x then result := result - 1;
  end else result := x;
end;

function nearrec(x: numeric): numeric;
begin
  result := nearint(recip(x));
  if isfin(result) and (result = 0)
    then result := nearint(x)
    else result := recip(result);
end;

function floor(x: numeric): ordinal;
begin
  result := trunc(x);
  if frac(x) < 0 then result := result - 1;
end;

function ceil(x: numeric): ordinal;
begin
  result := -floor(-x);
end;

function round(x: numeric): ordinal;
begin
  result := floor(x + 1/2);
end;

function modulo(x: numeric; y: numeric = 1): numeric;
begin
  result := x - y * floor(x/y);
end;

function wrap(x, l, r: numeric): numeric;
begin
  result := modulo(x - l, r - l) + l;
end;

function wrap(x, y: numeric): numeric;
begin
  result := wrap(x, y, -y);
end;

function wrap(x: numeric): numeric;
begin
  result := wrap(x, π);
end;

function deg(rad: numeric): numeric;
begin
  result := rad * 180/π;
end;

function rad(deg: numeric): numeric;
begin
  result := deg * π/180;
end;


function sqr(x: numeric): numeric;
begin
  result := x * x;
end;

function cub(x: numeric): numeric;
begin
  result := x * x * x;
end;

{$IFDEF WIRTH}

function sqrt(x: numeric): numeric;
var
  s: numeric;
begin
  if isnan(x) or (x < 0) then result := nan else
  if (x = 0) or (x = 1) or (x = inf) then result := x else
  begin
    s := x;  result := power(2, (ilog2(s) + 1) div 2);
    repeat
      s := result;  result := amean(s, x/s);
    until result = s;
  end;
end;

function cbrt(x: numeric): numeric;
var
  s: numeric;
begin
  if isnan(x) then result := nan else begin
    s := abs(x);
    if (s = 0) or (s = 1) or (s = inf) then result := x else
    begin
      result := sign(x, power(2, (ilog2(s) + 2) div 3));
      repeat
        s := result;  result := (2*s + x/(s * s))/3;
      until result = s;
    end;
  end;
end;

{$ELSE}


function sqrt(x: numeric): numeric;
asm
  FLD x; FSQRT; FWAIT
end;

function cbrt(x: numeric): numeric;
begin
  if x = 0 then result := 0 else
    result := sign(x, exp2(log2(abs(x))/3));
end;

{$ENDIF}

function hypot(x, y: numeric): numeric;
begin
  x := abs(x);  y := abs(y);
  if x > y then exchange(x, y);
  if x = 0
    then result := y
    else result := y * sqrt(1 + sqr(x/y));
end;

function compmod(x: numeric): numeric;
begin
  result := sqrt(1 - sqr(x));
end;

function TimeDilation(v: numeric): numeric;
const
  ϲ = 299792458;
begin
  result := compmod(v/ϲ);
end;

{$IFDEF WIRTH}

function arctan(x: numeric): numeric;
const
  maxiter = 64 * 2;
  halving = 1/256; // [0, 1)
var
  a, q, s, c, p: numeric;
  m, n: integer;
begin
  s := abs(x);
  if s = 0 then result :=  0  else
  if s = 1 then result := π/4 else
  if s = ∞ then result := π/2 else
  begin
    if s > 1 then begin
      s := 1/s;  x := -x;  p := π/2;
    end else p := 0;
    q := 1 + s * s;  m := 1;
    if q > 1 then begin
      if s > halving then begin
        c := (q + 1)/2;
        repeat
          a := c;  c := (a + q/a)/2;
        until c = a;
        s := s / (1 + c);
        q := 1 + s * s;  m := 2;
      end;
      n := 1;  a := s / q;  q := s * a;
      result := a;  c := 0;
      repeat
        a := a * q;
        n := n + 1;  a := a * n;
        n := n + 1;  a := a / n;
      until Kahan(result, s, c, a) or (n >= maxiter);
    end else result := s; // out of precision
    result := result * m - p;
  end;
  if x < 0 then result := -result;
end;

function arg(x, y: numeric): numeric;
begin
  if x = 0 then result := sign(y, π/2) else
  begin
    result := arctan(y/x);
    if x < 0 then result := result + psign(y, π);
  end;
end;

{$ELSE}

function arg(x, y: numeric): numeric;
asm
  FLD y; FLD x; FPATAN; FWAIT
end;

function arctan(x: numeric): numeric;
begin
  result := arg(1, x);
end;

{$ENDIF}

function atan2(y, x: numeric): numeric;
begin
  result := arg(x, y);
end;

function arctanloc(x, location: numeric): numeric;
begin
  result := arctan(x);
  result := result + π * round((location - result)/π);
end;

function arcsin(x: numeric): numeric;
begin
  if x < 0 then result := -arcsin(-x) else
  if x > 1 then result := Ø else
  if x = 0 then result := 0 else
  if x = 1 then result := π/2 else
  if x = 1/2 then result := π/6 else
  if x = sqrth then result := π/4 else
  if x = sqrt3/2 then result := π/3 else
  if x = (φ - 1)/2 then result := π/10 else
    result := 2 * arg(1 + compmod(x), x);
end;

function arccos(x: numeric): numeric;
begin
  result := π/2 - arcsin(x);
end;

function arccot(x: numeric): numeric;
begin
  result := arg(x, 1);
end;

function arccsc(x: numeric): numeric;
begin
  result := arcsin(recip(x));
end;

function arcsec(x: numeric): numeric;
begin
  result := π/2 - arccsc(x);
end;

function archav(x: numeric): numeric;
begin
  if x < 0 then result := Ø else
    result := 2 * arcsin(sqrt(x));
end;

function arccrd(x: numeric): numeric;
begin
  result := 2 * arcsin(x/2);
end;


{$IFDEF WIRTH}

function lniter(s, q: numeric): numeric; overload;
const
  maxiter = 64 * 2 + 1;
var
  n: integer;
  a, r, c: numeric;
begin
  n := 1;  a := s;  c := 0;  result := a;
  repeat
    n := n + 2;  a := a * q;
  until Kahan(result, r, c, a/n) or (n > maxiter);
  result := 2 * result;
end;

function lniter(s: numeric): numeric; overload;
begin
  result := lniter(s, s*s);
end;

function log2(x: numeric): numeric;
var
  l: ordinal;
begin
  if x < 0 then result := Ø else
  if x = 0 then result := -∞ else
  if x = ∞ then result := x else
  if x = ℮ then result := 1/ln2 else
  if x = 10 then result := lg10 else
  begin
    l := ilog2(x);
    if x = 1 then result := 0 else
    if x = sqrt2 then result := 1/2 else
      result := lniter((x - 1)/(x + 1)) / ln2;
    result := result + l;
  end;
end;

{$ELSE}

function log2(x: numeric): numeric;
asm
  FLD1; FLD x; FYL2X; FWAIT
end;

{$ENDIF}


function log(x, y: numeric): numeric;
begin
  if x = y then result := 1 else
  if y = 1 then result := 0 else
  begin
    result := log2(y);
    if (x <> 2) and (result <> 0) and (abs(result) <> ∞) then
      if x =  ℮ then result := result * ln2 else
      if x = 10 then result := result / lg10 else
        result := result / log2(x);
  end;
end;

function ln(x: numeric): numeric;
begin
  result := log(℮, x);
end;

function log10(x: numeric): numeric;
begin
  result := log(10, x);
end;

function lnone(x: numeric): numeric;
begin
  result := -ln(1 - x);
end;

function naplog(x: numeric; n: numeric = 10000000): numeric;
begin
  result := -n * ln(x/n);
end;

function napexp(x: numeric; n: numeric = 10000000): numeric;
begin
  result := n / exp(x/n);
end;

function mflog(x: numeric): numeric;
begin
  result := ln(x) + 8*ln2;
end;

function mfexp(x: numeric): numeric;
begin
  result := exp(x)/256;
end;

function midifreq(note: numeric): numeric; // middle C4 = 60
begin
  result := 440 * exp2((note - 69)/12);
end;

function midinote(freq: numeric): numeric;
begin
  if freq <= 0 then result := -∞ else
    result := log2(freq/440) * 12 + 69;
end;

function pianofreq(key: numeric): numeric; // middle C4 = 40
begin
  result := midifreq(key + 20);
end;

function pianokey(freq: numeric): numeric;
begin
  if freq <= 0 then result := -∞ else
    result := midinote(freq) - 20;
end;


function power(x: numeric; n: ordinal): numeric;
var
  i: ordinal;
begin
  if (n = 0) then result := 1 else
  if (x = 0) or (x = 1) then result := x else
  begin
    result := 1;  i := abs(n);
    while i > 0 do begin
      if odd(i) then result := result * x;
      i := i div 2;
      if i > 0 then x := sqr(x);
    end;
  end;
  if n < 0 then result := recip(result);
end;


{$IFDEF WIRTH}

procedure taylorexp(x: numeric; var sin, cos, sinh, cosh, exp: numeric; maxiter: integer = 32 * 4);
var
  s, t, c: array [0..3] of numeric;
  f: array [0..3] of boolean;
  i, n: integer;
  a: numeric;
begin
  for i := 0 to 3 do begin
    s[i] := 0; c[i] := 0; t[i] := ∞; f[i] := false;
  end;
  n := 0; a := 1; s[n] := a;
  if x <> 0 then
  while (n < maxiter) and not (f[0] and f[1] and f[2] and f[3]) do
  begin
    n := n + 1;  a := a * x/n;
    i := n mod 4;  f[i] := Kahan(s[i], t[i], c[i], a);
  end;
  cos := s[0] - s[2];  cosh := s[0] + s[2];
  sin := s[1] - s[3];  sinh := s[1] + s[3];
  exp := s[0] + s[1] + s[2] + s[3];
end;

function exp2(x: numeric): numeric;
var
  f: numeric;
  n: ordinal;
begin
  if x = 0 then result := 1 else
  begin
    f := abs(x);  n := trunc(f);  f := frac(f);
    if f =  0  then result := 1 else
    if f = 1/2 then result := sqrt2 else
      taylorexp(f * ln2, f, f, f, f, result);
    result := result * power(2, n);
    if x < 0 then result := recip(result);
  end;
end;

{$ELSE}

function exp2(x: numeric): numeric;
asm
  FLD x             // x
  FLD ST(0)         // x, x
  FRNDINT           // x, trunc(x)
  FSUB ST(1), ST(0) // x - trunc(x) = frac(x), trunc(x)
  FXCH              // trunc(x), frac(x)
  F2XM1             // trunc(x), 2^frac(x) - 1
  FLD1              // trunc(x), 2^frac(x) - 1, 1
  FADD              // trunc(x), 2^frac(x)
  FSCALE            // trunc(x), 2^(frac(x) + trunc(x)) = 2^x
  FSTP ST(1)        // 2^x
  FWAIT
{
  FLD x             // x
  FLD1              // x, 1
  FLD ST(1)         // x, 1, x
  FPREM             // x, 1, x mod 1 = frac(x)
  FSUB ST(2), ST(0) // x - frac(x) = trunc(x), 1, frac(x)
  F2XM1             // trunc(x), 1, 2^frac(x) - 1
  FADD              // trunc(x), 2^frac(x)
  FSCALE            // trunc(x), 2^(frac(x) + trunc(x)) = 2^x
  FSTP ST(1)        // 2^x
  FWAIT
}
end;

{$ENDIF}

function power(x, y: numeric): numeric;
begin
  if (y = 0) or  (x = 1) then result := 1 else
  if (x = 0) and (y > 0) then result := 0 else
  if y =   1  then result := x else
  if y =  -1  then result := recip(x) else
  if y =   2  then result := sqr(x) else
  if y =  -2  then result := recip(sqr(x)) else
  if y =  1/2 then result := sqrt(x) else
  if y = -1/2 then result := recip(sqrt(x)) else
  if frac(y) = 0 then result := power(x, trunc(y)) else
  if x =   2  then result := exp2(y) else
  if x =   ℮  then result := exp(y) else
  if x =  10  then result := exp10(y) else
    result := exp2(y * log2(x));
end;

function exp(x: numeric): numeric;
begin
  result := exp2(x / ln2);
end;

function exp10(x: numeric): numeric;
begin
  result := power(10, trunc(x));  x := frac(x);
  if x <> 0 then result := result * exp2(x * lg10);
end;

function logit(x: numeric; a: numeric = 1): numeric;
begin
  if (x <= 0) then result := -∞ else
  if (x >= 1) then result :=  ∞ else
    result := -ln(power(x, -recip(a)) - 1);
end;

function logistic(x: numeric; a: numeric = 1): numeric;
begin
  if x = -∞ then result := 0 else
  if x =  ∞ then result := 1 else
    result := power(1 + exp(-x), -a);
end;

{$IFDEF WIRTH}

procedure sincos(x: numeric; var sin, cos: numeric);
var
  a: numeric;
begin
  if abs(x) > π then x := wrap(x);
  a := abs(x);
  if a = 0   then begin sin := 0;       cos :=  1 end else
  if a = π/2 then begin sin := sign(x); cos :=  0 end else
  if a = π   then begin sin := 0;       cos := -1 end else
    taylorexp(x, sin, cos, a, a, a);
end;

{$ELSE}

procedure sincos(x: numeric; var sin, cos: numeric);
asm
  FLD x; FSINCOS
  FSTP tbyte ptr [edx]
  FSTP tbyte ptr [eax]
  FWAIT
end;

{$ENDIF}


function sin(x: numeric): numeric;
begin
  sincos(x, result, x);
end;

function cos(x: numeric): numeric;
begin
  sincos(x, x, result);
end;

function tan(x: numeric): numeric;
var
  sin, cos: numeric;
begin
  sincos(x, sin, cos);  result := sin/cos;
end;

function cot(x: numeric): numeric;
var
  sin, cos: numeric;
begin
  sincos(x, sin, cos);  result := cos/sin;
end;

function sec(x: numeric): numeric;
begin
  result := recip(cos(x));
end;

function csc(x: numeric): numeric;
begin
  result := recip(sin(x));
end;

function hav(x: numeric): numeric;
begin
  result := sqr(sin(x/2));
end;

function crd(x: numeric): numeric;
begin
  result := 2 * sin(x/2);
end;

function sinc(x: numeric): numeric;
begin
  if x = 0
    then result := 1
    else result := sin(x)/x;
end;


procedure sindcosd(x: numeric; var sin, cos: numeric);
begin
  sincos(rad(x), sin, cos);
end;

function sind(x: numeric): numeric;
begin
  result := sin(rad(x));
end;

function cosd(x: numeric): numeric;
begin
  result := cos(rad(x));
end;

function tand(x: numeric): numeric;
begin
  result := tan(rad(x));
end;

function cotd(x: numeric): numeric;
begin
  result := cot(rad(x));
end;

function argd(x, y: numeric): numeric;
begin
  result := deg(arg(x, y));
end;

function arcsind(x: numeric): numeric;
begin
  result := deg(arcsin(x));
end;

function arccosd(x: numeric): numeric;
begin
  result := deg(arccos(x));
end;

function arctand(x: numeric): numeric;
begin
  result := deg(arctan(x));
end;

function arccotd(x: numeric): numeric;
begin
  result := deg(arccot(x));
end;


procedure sinhcosh(x: numeric; var sinh, cosh: numeric);
var
  p, q: numeric;
begin
  if x = 0 then begin
    sinh := 0; cosh := 1;
  end else begin
    if x < 0 then begin
      q := exp(-x); p := recip(q);
    end else begin
      p := exp(x); q := recip(p);
    end;
    cosh := amean(p, q);
    sinh := p - cosh;
  end;
end;

function sinh(x: numeric): numeric;
begin
  sinhcosh(x, result, x);
end;

function cosh(x: numeric): numeric;
begin
  sinhcosh(x, x, result);
end;

function tanh(x: numeric): numeric;
var
  sinh, cosh: numeric;
begin
  if abs(x) = ∞ then result := sign(x) else begin
    sinhcosh(x, sinh, cosh);  result := sinh/cosh;
  end;
end;

function coth(x: numeric): numeric;
begin
  result := recip(tanh(x));
end;

function arctanh(x: numeric): numeric;
begin
  result := abs(x);
  if result > 1 then result := Ø else
  if result = 1 then result := sign(x, ∞) else
  if x <> 0 then result := ln((1 + x)/(1 - x))/2;
end;

function arccoth(x: numeric): numeric;
begin
  if x = 0 then result := zsign(x, ∞) else
  if abs(x) = ∞ then result := sign(x, 0) else
    result := ln((x + 1)/(x - 1))/2;
end;

function arcsinh(x: numeric): numeric;
begin
  if x = 0 then result := 0 else
    result := ln(x + sqrt(sqr(x) + 1));
end;

function arccosh(x: numeric): numeric;
begin
  if x < 1 then result := Ø else
  if x = 1 then result := 0 else
    result := ln(x + sqrt(sqr(x) - 1));
end;

function hadd(x, y: numeric): numeric; // 1 / (1/x + 1/y)
begin
  if x = y then result := x/2 else
  if abs(x) = ∞ then result := sign(x, y) else
  if abs(y) = ∞ then result := sign(y, x) else
  if (x = 0) or (y = 0) then result := 0 else
  begin
    result := x + y;
    if result = 0
      then result := sign(x * y, ∞)
      else result := (x * y)/result;
  end;
end;

function hsub(x, y: numeric): numeric;
begin
  result := hadd(x, -y);
end;

function hmean(x, y: numeric): numeric;
begin
  result := 2 * hadd(x, y);
end;

function amean(x, y: numeric; t: numeric = 1/2): numeric;
begin
  if x = y then result := x else
    result := x + t * (y - x);
end;

function gmean(x, y: numeric): numeric;
begin
  if x = y then result := x else
    result := sqrt(x * y);
end;

function agm(x, y: numeric): numeric; // arithmetic-geometric mean
const
  maxiter = 64;
var
  i: integer;
  z: numeric;
begin
  if x = y then result := x else begin
    result := y;  y := 0;  i := 0;
    while (result <> x) and (result <> y) and (i < maxiter) do begin
      if x = y then result := x;
      y := result;  z := x;  i := i + 1;
      result := amean(x, y); x := gmean(x, y);
      if x = z then result := x;
    end;
  end;
end;

function EllipticK(x: numeric): numeric;
begin
  if x > 1 then result := Ø else
  if x = 1 then result := ∞ else
  if x = -∞ then result := 0 else
    result := π/2 / agm(1, sqrt(1 - x));
end;

function PendulumAmplitudePeriod(length, inclination, gravity: numeric): numeric;
begin
  result := 4 * sqrt(abs(length/gravity)) * EllipticK(hav(inclination));
end;

function Gudermannian(x: numeric): numeric;
begin
  result := 2 * arctan(tanh(x/2));
end;

function AntiCollision(n, k:integer): numeric;
begin
  result := 1;
  while k > 1 do begin
    k := k - 1;
    result := result * (n - k)/n;
  end;
end;

function Gaussian(x, a, b, c, d, p: numeric): numeric;
begin
  result := a * exp(-power((x - b)/c, p)/2) + d;
end;


function mvhgdprob(const play, hits: array of integer): numeric; overload;
var
  balls, draws, i: integer;
begin
  if length(play) = length(hits) then begin
    result := 1;  balls := 0;  draws := 0;  i := low(play);
    while (result <> 0) and (i <= high(play)) do
    if (play[i] < hits[i]) or (hits[i] < 0) then
      result := 0
    else begin
      result := result * binom(play[i], hits[i]);
      inc(balls, play[i]);  inc(draws, hits[i]);
      i := i + 1;
    end;
    if result <> 0 then result := result / binom(balls, draws);
  end else result := Ø;
end;

function hgdprob(balls, draws, play, hits: integer): numeric; overload;
begin
  result := mvhgdprob([play, balls - play], [hits, draws - hits]);
end;


function KenoProb(play, hits: integer): numeric; overload;
begin
  result := hgdprob(80, 20, play, hits);
end;


class operator dot.implicit(d: dot): complex;
begin
  result := xiy(d.x, d.y);
end;

class operator dot.implicit(z: complex): dot;
begin
  result.x := z.Χ;  result.y := z.Υ;
end;


function xiy(x: numeric; y: numeric = 0): complex;
begin
  result.x := x;
  result.y := y;
end;

function cis(θ: numeric): complex;
begin
  sincos(θ, result.y, result.x);
end;

function cis(θ: numeric; var x, y: numeric): complex;
begin
  result := cis(θ);  x := result.x;  y := result.y;
end;

function taxicab(z: complex): numeric;
begin
  result := abs(z.x) + abs(z.y);
end;

function sqrabs(z: complex): numeric;
begin
  result := sqr(z.x) + sqr(z.y);
end;

function abs(z: complex): numeric;
begin
  if (z.x = 0) or (z.y = 0)
    then result := taxicab(z)
    else result := sqrt(sqrabs(z));
end;

function sign(z: complex): complex;
begin
  if z = 0 then result := z else result := z/abs(z);
end;

function sign(x: numeric; z: complex): complex;
begin
  if x < 0 then result := -z else
  if x > 0 then result :=  z else
    result := 0;
end;

function arg(z: complex): numeric;
begin
  result := arg(z.x, z.y);
end;

function cisd(θ: numeric): complex;
begin
  result := cis(rad(θ));
end;

function argd(z: complex): numeric;
begin
  result := deg(arg(z));
end;

function eis(z: complex): complex;
begin
  z.y := z.y/2;  result := xiy(z.x - z.y, z.y * sqrt3);
end;

function complex.getabs: numeric;
begin
  result := abs(self);
end;

function complex.getarg: numeric;
begin
  result := arg(self);
end;

function complex.getreal: boolean;
begin
  result := y = 0;
end;

procedure complex.setreal(f: boolean = true);
begin
  if f then y := 0;
end;

function complex.getimag: boolean;
begin
  result := x = 0;
end;

procedure complex.setimag(f: boolean = true);
begin
  if f then x := 0;
end;

function complex.getintX: ordinal;
begin
  result := round(x);
end;

function complex.getintY: ordinal;
begin
  result := round(y);
end;


class operator complex.implicit(x: numeric): complex;
begin
  result := xiy(x);
end;

class operator complex.equal(u, v: complex): boolean;
begin
  result := (u.x = v.x) and (u.y = v.y);
end;

class operator complex.notequal(u, v: complex): boolean;
begin
  result := not (u = v);
end;

class operator complex.positive(z: complex): complex;
begin
  result := xiy(z.x, z.y);
end;

class operator complex.negative(z: complex): complex;
begin
  result := xiy(-z.x, -z.y);
end;

class operator complex.logicalnot(z: complex): complex;
begin
  result := xiy(z.x, -z.y);
end;

class operator complex.add(u, v: complex): complex;
begin
  result := xiy(u.x + v.x, u.y + v.y);
end;

class operator complex.add(z: complex; x: numeric): complex;
begin
  result := xiy(z.x + x, z.y);
end;

class operator complex.add(x: numeric; z: complex): complex;
begin
  result := z + x;
end;

class operator complex.subtract(u, v: complex): complex;
begin
  result := -v + u;
end;

class operator complex.subtract(z: complex; x: numeric): complex;
begin
  result := -x + z;
end;

class operator complex.subtract(x: numeric; z: complex): complex;
begin
  result := -z + x;
end;

class operator complex.multiply(u, v: complex): complex;
begin
  result := xiy(u.x * v.x - u.y * v.y, u.x * v.y + u.y * v.x);
end;

class operator complex.multiply(z: complex; x: numeric): complex;
begin
  result := xiy(z.x * x, z.y * x);
end;

class operator complex.multiply(x: numeric; z: complex): complex;
begin
  result := z * x;
end;

class operator complex.divide(z: complex; x: numeric): complex;
begin
  if x = 0
    then result := complexinf
    else result := xiy(z.x/x, z.y/x);
end;

class operator complex.divide(x: numeric; z: complex): complex;
begin
  result := x * (not z) / sqrabs(z);
end;

class operator complex.divide(u, v: complex): complex;
begin
  result := u * (not v) / sqrabs(v);
end;


function dotcross(u, v: complex): complex;
begin
  result := (not u) * v;
end;

function dotprod(u, v: complex): numeric;
begin
  result := dotcross(u, v).x;
end;

function crossprod(u, v: complex): numeric;
begin
  result := dotcross(u, v).y;
end;

function intersection(u1, u2, v1, v2: complex): complex;
begin
  u2 := u2 - u1;  v2 := v2 - v1;
  result := (crossprod(u2, u1) * v2 - crossprod(v2, v1) * u2)/crossprod(u2, v2);
end;


function floor(z: complex): complex;
begin
  result := xiy(floor(z.x), floor(z.y));
end;

function ceil(z: complex): complex;
begin
  result := xiy(ceil(z.x), ceil(z.y));
end;

function round(z: complex): complex;
begin
  result := xiy(round(z.x), round(z.y));
end;

function modulo(u, v: complex): complex;
begin
  result := u - v * floor(u/v);
end;


function sqr(z: complex): complex;
begin
  result := xiy((z.x + z.y) * (z.x - z.y), 2 * z.x * z.y);
end;

function sqrt(z: complex): complex;
begin
  if z.isreal then
    if z.x < 0
      then result := xiy(0, sqrt(-z.x))
      else result := sqrt(z.x)
  else
  if z.isimag then begin
    result.x := sqrt(abs(z.y)/2);
    result.y := sign(z.y, result.x);
  end else begin
    result.x := sqrt(amean(z.ρ, z.x));
    result.y := z.y / result.x / 2;
  end;
end;

function compmod(z: complex): complex;
begin
  result := sqrt(1 - sqr(z));
end;


function exp(z: complex): complex;
begin
  if z = complexinf then result := 0 else begin
    result := exp(z.x);
    if not z.isreal then result := result * cis(z.y);
  end;
end;

function ln(z: complex): complex;
begin
  if z = 0 then result := complexinf else
    result := xiy(ln(z.ρ), z.θ);
end;

function power(z: complex; n: ordinal): complex;
var
  i: ordinal;
begin
  if (n = 0) then result := 1 else
  if (z = 0) or (z = 1) then result := z else
  begin
    result := 1; i := abs(n);
    while i > 0 do begin
      if odd(i) then result := result * z;
      i := i div 2;
      if i > 0 then z := sqr(z);
    end;
  end;
  if n < 0 then result := 1/result;
end;

function power(z: complex; x: numeric): complex;
begin
  if z.isreal and (z.x >= 0) then result := power(z.x, x) else
  if x =    0 then result := 1 else
  if x =    1 then result := z else
  if x =   -1 then result := 1/z else
  if x =    2 then result := sqr(z) else
  if x =   -2 then result := 1/sqr(z) else
  if x =  1/2 then result := sqrt(z) else
  if x = -1/2 then result := 1/sqrt(z) else
    result := exp(x * ln(z));
end;

function power(u, v: complex): complex;
begin
  if v.isreal then result := power(u, v.x) else
    result := exp(v * ln(u));
end;

function power(x: numeric; z: complex): complex;
begin
  if z.isreal then result := power(x, z.x) else
    result := power(xiy(x), z);
end;

function log2(z: complex): complex;
begin
  if z.isreal and (z.x > 0) then result := log2(z.x) else
    result := ln(z) / ln2;
end;

function exp2(z: complex): complex;
begin
  if z.isreal then result := exp2(z.x) else
    result := exp(z * ln2);
end;

function root(z: complex; n: ordinal): complex;
begin
  result := power(z.ρ, 1/n) * cis(z.θ/n);
end;

function amean(u, v: complex; t: numeric = 1/2): complex;
begin
  if u = v then result := u else
    result := u + t * (v - u);
end;

function curve(u, p, v: complex; t: numeric = 1/2): complex;
begin
  result := amean(amean(u, p, t), amean(p, v, t), t);
end;

function bezier(u, p, q, v: complex; t: numeric = 1/2): complex;
begin
  result := amean(curve(u, p, q, t), curve(p, q, v, t), t);
end;

function gmean(u, v: complex): complex; overload;
begin
  if u = v then result := u else begin
    result := u * v;
    if result <> 0 then
    begin
      result := sqrt(result);
      result := sign(dotprod(result, amean(u, v)), result);
    end;
  end;
end;

function agm(u, v: complex): complex;
const
  maxiter = 64;
var
  i: integer;
  z: complex;
begin
  if u = v then result := u else begin
    result := v;  v := 0;  i := 0;
    while (result <> u) and (result <> v) and (i < maxiter) do begin
      if u = v then result := u;
      v := result;  z := u;  i := i + 1;
      result := amean(u, v); u := gmean(u, v);
      if u = z then result := u;
    end;
  end;
end;


procedure sinhcosh(z: complex; var sinh, cosh: complex);
var
  p, q: complex;
begin
  if z.isreal then begin
    sinhcosh(z.x, sinh.x, cosh.x);
    sinh.y := 0;  cosh.y := 0;
  end else begin
    if z.y < 0 then begin
      q := exp(-z); p := 1/q;
    end else begin
      p := exp(z); q := 1/p;
    end;
    cosh := amean(p, q);
    sinh := p - cosh;
  end;
end;

function sinh(z: complex): complex;
begin
  sinhcosh(z, result, z);
end;

function cosh(z: complex): complex;
begin
  sinhcosh(z, z, result);
end;

function tanh(z: complex): complex;
var
  sinh, cosh: complex;
begin
  sinhcosh(z, sinh, cosh);  result := sinh/cosh;
end;

function coth(z: complex): complex;
var
  sinh, cosh: complex;
begin
  sinhcosh(z, sinh, cosh);  result := cosh/sinh;
end;

procedure sincos(z: complex; var sin, cos: complex);
begin
  if z.isreal then begin
    sincos(z.x, sin.x, cos.x);
    sin.y := 0;  cos.y := 0;
  end else begin
    sinhcosh(-z * ί, sin, cos);
    sin := sin * ί;
  end;
end;

function sin(z: complex): complex;
begin
  sincos(z, result, z);
end;

function cos(z: complex): complex;
begin
  sincos(z, z, result);
end;

function tan(z: complex): complex;
var
  sin, cos: complex;
begin
  sincos(z, sin, cos);  result := sin/cos;
end;

function cot(z: complex): complex;
var
  sin, cos: complex;
begin
  sincos(z, sin, cos);  result := cos/sin;
end;

function sec(z: complex): complex;
begin
  result := 1/cos(z);
end;

function csc(z: complex): complex;
begin
  result := 1/sin(z);
end;

function hav(z: complex): complex;
begin
  result := sqr(sin(z/2));
end;

function crd(z: complex): complex;
begin
  result := 2 * sin(z/2);
end;

function sinc(z: complex): complex;
begin
  if z = 0
    then result := 1
    else result := sin(z)/z;
end;


function arcsinh(z: complex): complex;
begin
  result := ln(z + sqrt(sqr(z) + 1));
end;

function arccosh(z: complex): complex;
begin
  result := ln(z + gmean(z + 1, z - 1));
end;

function arctanh(z: complex): complex;
begin
  if z.isreal and (abs(z) < 1)
    then result := arctanh(z.x)
    else result := ln((1 + z)/(1 - z)) / 2;
end;

function arccoth(z: complex): complex;
begin
  if z = complexinf
    then result := 0
    else result := ln((z + 1)/(z - 1)) / 2;
end;


function arcsin(z: complex): complex;
begin
  if z.isreal and (abs(z.x) <= 1)
    then result := arcsin(z.x)
    else result := -ί * arcsinh(z * ί);
end;

function arccos(z: complex): complex;
begin
  result := π/2 - arcsin(z);
end;

function arctan(z: complex): complex;
begin
  if z.isreal
    then result := arctan(z.x)
    else result := -ί * arctanh(z * ί);
end;

function arccot(z: complex): complex;
begin
  if z.isreal
    then result := arccot(z.x)
    else result := π/2 - arctan(z);
end;

function arccsc(z: complex): complex;
begin
  result := arcsin(1/z);
end;

function arcsec(z: complex): complex;
begin
  result := π/2 - arccsc(z);
end;

function archav(z: complex): complex;
begin
  result := 2 * arcsin(sqrt(z));
end;

function arccrd(z: complex): complex;
begin
  result := 2 * arcsin(z/2);
end;

function arctanloc(z, location: complex): complex;
begin
  result := arctan(z);
  result := result + π * round((location - result)/π);
end;


function hms(h: numeric; m: numeric = 0; s: numeric = 0): numeric;
begin
  result := psign(h, abs(h) + (m + s / 60) / 60);
end;

function dms(d: numeric; m: numeric = 0; s: numeric = 0): numeric;
begin
  result := rad(hms(d, m, s));
end;

function geopos(lat, lon: numeric): complex;
begin
  lat := wrap(lat, 180);
  lon := wrap(lon, 180);
  if abs(lat) > 90 then begin
    lat := sign(lat, 180 - abs(lat));
    lon := lon - psign(lon, 180);
  end;
  result := xiy(rad(lat), rad(lon));
end;

function geodist(a, b: complex; r: numeric): numeric;
begin
  result := r * archav(hav(a.x - b.x) + hav(a.y - b.y) * cos(a.x) * cos(b.x));
end;

function mandelbrot(c, p: complex; m: integer = 1000): integer;
var
  z, w: complex;
begin
  // result: < 0: in set; = 0: probably in set; > 0: not in set; = m: out of bounds.
  result := m;  z := c;
  while (result > 0) and (sqrabs(z) <= 4) do begin
    w := z;  z := power(z, p) + c;
    if z = w
      then result := -result
      else result := result - 1;
  end;
end;

function mandelbrot(c: complex; x: numeric = 2; m: integer = 1000): integer;
begin
  result := mandelbrot(c, xiy(x), m);
end;

function superellipse(radius: complex; angle, shape, symmetry, u, v: numeric): complex; overload;
var
  z, w: complex;
  r: numeric;
begin
  if radius = 0 then result := 0 else
  begin
    w := cis(angle);  symmetry := abs(symmetry);
    if radius.y = 0 then result := w.x * radius else
    if radius.x = 0 then result := w.y * radius else
    begin
      if symmetry = 4
        then z := w
        else z := cis(angle * symmetry/4);
      z := xiy(abs(z.x/radius.x), abs(z.y/radius.y));
      z := xiy(power(z.x, u), power(z.y, v));
      r := power(taxicab(z), -recip(shape));
      result := r * w;
    end;
  end;
end;

function superellipse(radius: complex; angle, shape, symmetry: numeric): complex; overload;
begin
  result := superellipse(radius, angle, shape, symmetry, shape, shape);
end;

function superellipse(radius: numeric; angle, shape, symmetry: numeric): complex; overload;
begin
  result := superellipse(xiy(radius, radius), angle, shape, symmetry);
end;

function PendulumAmplitudePeriod(body: complex; gravity: numeric = 9.80665): numeric;
begin
  result := PendulumAmplitudePeriod(body.ρ, body.θ, gravity);
end;


function poly(x: numeric; const a: array of numeric): numeric;
var
  i: integer;
begin
  result := 0;
  for i := high(a) downto low(a) do result := result * x + a[i];
end;

function poly(z: complex; const a: array of numeric): complex;
var
  i: integer;
begin
  result := 0;
  for i := high(a) downto low(a) do result := result * z + a[i];
end;

function poly(z: complex; const a: array of complex): complex;
var
  i: integer;
begin
  result := 0;
  for i := high(a) downto low(a) do result := result * z + a[i];
end;

function TropicalYear(year: numeric): numeric;
begin
  result := poly(year - 2000, [365.242189669781, 6.16187e-6, -6.44e-10]);
end;


// Taylor Serie - Gamma Reciprocal
// N[CoefficientList[Series[1/Gamma[z], {z, 0, 30}], z], 24]
const
  tsgr: array [0..31] of numeric = (0, 1,
     0.57721566490153286060651209008240243104215933593992, // γ (Euler–Mascheroni constant)
    -0.6558780715202538810770195151453904812798, // (γ² - π²/6)/2
    -0.0420026350340952355290039348754298187114,
     0.1665386113822914895017007951021052357178,
    -0.0421977345555443367482083012891873913017,
    -9.62197152787697356211492e-03,
     7.21894324666309954239501e-03,
    -1.16516759185906511211397e-03,
    -2.15241674114950972815730e-04,
     1.28050282388116186153199e-04,
    -2.01348547807882386556894e-05,
    -1.25049348214267065734536e-06,
     1.13302723198169588237413e-06,
    -2.05633841697760710345015e-07,
     6.11609510448141581786250e-09,
     5.00200764446922293005567e-09,
    -1.18127457048702014458813e-09,
     1.04342671169110051049154e-10,
     7.78226343990507125404994e-12,
    -3.69680561864220570818782e-12,
     5.10037028745447597901548e-13,
    -2.05832605356650678322243e-14,
    -5.34812253942301798237002e-15,
     1.22677862823826079015889e-15,
    -1.18125930169745876951376e-16,
     1.18669225475160033257978e-18,
     1.41238065531803178155580e-18,
    -2.29874568443537020659248e-19,
     1.71440632192733743338396e-20,
     1.33735173049369311486478e-22);

function Gamma(x: numeric): numeric;
begin
  if abs(x) > 32 then begin
    result := exp2(x - 1);  x := x/2;
    result := Gamma(x) * Gamma(x + 1/2) * result / sqrtπ;
  end else begin
    result := 1;
    while x >  1 do begin x := x - 1; result := result * x end;
    while x <= 0 do begin result := result / x; x := x + 1 end;
    if x <> 1 then
      if x = 1/2
        then result := result * sqrtπ
        else result := result / poly(x, tsgr);
  end;
end;

function fact(x: numeric): numeric;
begin
  if (x < 0) and (frac(x) = 0) then begin
    result := 1/Gamma(-x);
    if not odd(trunc(x)) then result := -result;
  end else result := Gamma(x + 1);
end;

function fact2(x: numeric): numeric;
begin
  result := x/2;
  result := fact(result) * exp2(result);
  if frac(x) = 0 then begin
    if odd(trunc(x)) then result := result * sqrt2/sqrtπ;
    if x < 0
      then result := nearrec(result)
      else result := nearint(result);
  end else result := result * power(π/2, (cos(π*x) - 1)/4)
end;

function binom(n, k: numeric): numeric;
var
  l: numeric;
begin
  l := n - k;
  if (k = 0) or (l = 0) then result := 1 else
  if (k = 1) or (l = 1) then result := n else
  begin
    result := fact(n);
    if k = l
      then result := result / sqr(fact(k))
      else result := result / (fact(k) * fact(l));
  end;
end;

function Gamma(z: complex): complex;
begin
  if z.isreal then result := Gamma(z.x) else
  if (abs(z.y) > 1) or (abs(z.x) > 32) then begin
    result := exp2(z - 1);  z := z/2;
    result := Gamma(z) * Gamma(z + 1/2) * result / sqrtπ;
  end else begin
    result := 1;
    while z.x >  1 do begin z := z - 1; result := result * z end;
    while z.x <= 0 do begin result := result / z; z := z + 1 end;
    if z <> 1 then
      if z = 1/2
        then result := result * sqrtπ
        else result := result / poly(z, tsgr);
  end;
end;

function fact(z: complex): complex;
begin
  if z.isreal
    then result := fact(z.x)
    else result := Gamma(z + 1);
end;

function fact2(z: complex): complex;
begin
  if z.isreal then result := fact2(z.x) else
    result := fact(z/2) * exp2(z/2) * power(π/2, (cos(π*z) - 1)/4);
end;

function binom(n, k: complex): complex;
var
  l: complex;
begin
  l := n - k;
  if (k = 0) or (l = 0) then result := 1 else
  if (k = 1) or (l = 1) then result := n else
  begin
    result := fact(n);
    if k = l
      then result := result / sqr(fact(k))
      else result := result / (fact(k) * fact(l));
  end;
end;

function multinom(const x: array of numeric): numeric; overload;
var
  s: numeric;
  i: integer;
begin
  s := 0;  result := 1;
  for i := low(x) to high(x) do begin
    s := s + x[i];  result := result * binom(s, x[i]);
  end;
end;

function multinom(const z: array of complex): complex; overload;
var
  s: complex;
  i: integer;
begin
  s := 0;  result := 1;
  for i := low(z) to high(z) do begin
    s := s + z[i];  result := result * binom(s, z[i]);
  end;
end;

function fibnumber(n: ordinal): numeric;
begin
  result := nearint(power(φ, abs(n))/sqrt5);
  if (n < 0) and odd(n) then result := -result;
end;

function fibnumber(x: numeric): complex;
begin
  if frac(x) = 0 then result := fibnumber(trunc(x)) else
    result := (power(φ, x) - power(xiy(-φ), -x))/sqrt5;
end;

procedure exchange(var a, b: ordinal); overload;
var
  c: ordinal;
begin
  c := a;  a := b;  b := c;
end;

procedure exchange(var a, b: numeric); overload;
var
  c: numeric;
begin
  c := a;  a := b;  b := c;
end;

procedure exchange(var a, b: complex); overload;
var
  c: complex;
begin
  c := a;  a := b;  b := c;
end;

function Beta(u, v: ordinal): numeric; overload;
var
  w: ordinal;
begin
  if u > v then exchange(u, v);  w := u + v;
  if w = 0 then begin
    if u = v then result := ∞ else begin
      result := 1/v;
      if odd(v) then result := -result;
    end
  end else
  if (u = 0) or (v = 0) then result := ∞ else
  if u < 0 then begin
    if v < 0 then result := ∞ else
    if (w > u) and (w < 0)
      then result := fact(u-1) * fact(v-1) / fact(w-1)
      else result := ∞;
  end else begin
    if u = v
      then result := sqr(fact(u-1))/fact(w-1)
      else result := fact(u-1) * fact(v-1) / fact(w-1)
  end;
end;

function Beta(u, v: numeric): numeric; overload;
var
  w: numeric;
begin
  if (frac(u) = 0) and (frac(v) = 0) then
    result := Beta(trunc(u), trunc(v))
  else begin
    if u > v then exchange(u, v);  w := u + v;
    if (u = 0) or (v = 0) then result := ∞ else
    if (frac(u) = 0) and (v < 0) then result := ∞ else
    if (frac(v) = 0) and (u < 0) then result := ∞ else
    if u = 1 then result := 1/v else
    if v = 1 then result := 1/u else
    if w = 0 then result := 0 else
    if w = 1 then
      if u = v
        then result := π
        else result := π / sin(π * v)
    else
    if u = v
      then result := sqr(fact(u-1))/fact(w-1)
      else result := fact(u-1) * fact(v-1) / fact(w-1)
  end;
end;

function Beta(u, v: complex): complex; overload;
var
  w: complex;
begin
  if u.isreal and v.isreal then
    result := Beta(u.x, v.x)
  else begin
    if u.x > v.x then exchange(u, v);  w := u + v;
    if (u = 0) or (v = 0) then result := complexinf else
    if u = 1 then result := 1/v else
    if v = 1 then result := 1/u else
    if w = 0 then result := 0 else
    if w = 1 then result := π / sin(π * v) else
    if u = v
      then result := sqr(fact(u-1))/fact(w-1)
      else result := fact(u-1) * fact(v-1) / fact(w-1)
  end;
end;

function fallfact(x, n: numeric): numeric; overload;
begin
  result := fact(x) / fact(x - n);
end;

function risefact(x, n: numeric): numeric; overload;
begin
  x := x - 1;  result := fact(x + n) / fact(x);
end;

function fallfact(z, n: complex): complex; overload;
begin
  result := fact(z) / fact(z - n);
end;

function risefact(z, n: complex): complex; overload;
begin
  z := z - 1;  result := fact(z + n) / fact(z);
end;


function BallVolume(dimension: numeric): numeric;
begin
  result := dimension/2;
  result := power(π, result) / fact(result);
end;

function SphereArea(dimension: numeric): numeric;
begin
  result := dimension * BallVolume(dimension);
end;

function Lissajous(t, a, b, u, v, d: numeric): complex;
begin
  result := xiy(a * sin(u * t + d), b * sin(v * t));
end;

function Lissajous(t, a, b, u, v: numeric): complex;
begin
  result := Lissajous(t, a, b, u, v, π/2);
end;

function Lissajous(t, u, v: numeric): complex;
begin
  result := Lissajous(t, 1, 1, u, v);
end;


function frc(p: ordinal; q: ordinal = 1): rational;
var
  g: ordinal;
begin
  if q = 0 then p := sign(p) else
  begin
    g := sign(q, abs(gcd(p, q)));
    p := p div g;  q := q div g;
  end;
  result.p := p;  result.q := q;
end;

function rational.getx: numeric;
begin
  if q = 0 then
    if p = 0
      then result := Ø
      else result := sign(p, ∞)
  else result := p/q;
end;

procedure rational.setx(r: numeric);
begin
  fractional(r, p, q);
end;

class operator rational.implicit(x: numeric): rational;
begin
  if isnan(x) then result := frc(0, 0) else
  if isinf(x) then result := frc(trunc(sign(x)), 0) else
    result.x := x;
end;

class operator rational.implicit(n: ordinal): rational;
begin
  result := frc(n);
end;

class operator rational.implicit(f: rational): numeric;
begin
  result := f.x;
end;

class operator rational.implicit(f: rational): complex;
begin
  result := f.x;
end;

class operator rational.equal(f, g: rational): boolean;
begin
  result := f.x = g.x;
end;

class operator rational.notequal(f, g: rational): boolean;
begin
  result := not (f = g);
end;

class operator rational.lessthan(f, g: rational): boolean;
begin
  result := f.x < g.x;
end;

class operator rational.greaterthan(f, g: rational): boolean;
begin
  result := f.x > g.x;
end;

class operator rational.lessthanorequal(f, g: rational): boolean;
begin
  result := not (f > g);
end;

class operator rational.greaterthanorequal(f, g: rational): boolean;
begin
  result := not (f < g);
end;

class operator rational.positive(f: rational): rational;
begin
  result := frc(f.p, f.q);
end;

class operator rational.negative(f: rational): rational;
begin
  result := frc(-f.p, f.q);
end;

class operator rational.logicalnot(f: rational): rational;
begin
  result := frc(f.q, f.p);
end;

class operator rational.multiply(f, g: rational): rational;
begin
  f := +f;  g := +g; // ensure reduced
  result.q := f.q;  f := frc(f.p, g.q);  g := frc(g.p, result.q);
  result := frc(f.p * g.p, f.q * g.q);
end;

class operator rational.divide(f, g: rational): rational;
begin
  result := not g * f;
end;

class operator rational.intdivide(f, g: rational): ordinal;
begin
  result := trunc(f/g);
end;

class operator rational.modulus(f, g: rational): rational;
begin
  result := f - g * floor(f/g);
end;

class operator rational.add(f, g: rational): rational;
var
  d: ordinal;
begin
  if (f.q = 0) or (g.q = 0) then begin
    result := frc(sign(f.p) * sign(g.p), 0);
  end else begin
    f := +f;  g := +g; // ensure reduced
    d := lcm(f.q, g.q);
    f.p := (d div f.q) * f.p; f.q := d;
    g.p := (d div g.q) * g.p; g.q := d;
    result := frc(f.p + g.p, d);
  end;
end;

class operator rational.subtract(f, g: rational): rational;
begin
  result := -g + f;
end;

class operator rational.inc(f: rational): rational;
begin
  result := f + 1;
end;

class operator rational.dec(f: rational): rational;
begin
  result := f - 1;
end;

class operator rational.trunc(f: rational): ordinal;
begin
  result := trunc(f);
end;

function gcd(f, g: rational): rational;
begin
  result := f;
  while g <> 0 do begin
    f := g;
    g := result mod g;
    result := f;
  end;
end;

function lcm(f, g: rational): rational;
begin
  if (f = 0) and (g = 0)
    then result := 0
    else result := f / gcd(f, g) * g;
end;

function sqr(f: rational): rational;
begin
  result := f * f;
end;

function amean(f, g: rational): rational;
begin
  result := (f + g)/2;
end;

function hadd(f, g: rational): rational;
begin
  result := not (not f + not g);
end;

function hsub(f, g: rational): rational;
begin
  result := hadd(f, -g);
end;

function hmean(f, g: rational): rational;
begin
  if f = g then result := f else
  if (f = 0) or (g = 0) then result := 0 else
    result := 2 * hadd(f, g);
end;

function trunc(f: rational): ordinal;
begin
  result := f.p div f.q;
end;

function frac(f: rational): rational;
begin
  result := f - trunc(f);
end;

function floor(f: rational): ordinal;
begin
  result := trunc(f);
  if result > f then dec(result);
end;

function ceil(f: rational): ordinal;
begin
  result := -floor(-f);
end;

function round(f: rational): ordinal;
begin
  result := floor(f + frc(1, 2));
end;

function sign(f: rational): ordinal;
begin
  if f.q < 0 then result := -sign(f.p) else result := sign(f.p);
end;

function contfrac(const v: array of ordinal; n: integer = -1): rational;
var
  i: integer;
begin
  if n < 0 then n := length(v) else n := min(n, length(v));
  if n = 0 then result := 0 else begin
    result.p := 1;  result.q := 0;
    for i := n - 1 + low(v) downto low(v) do
      result := frc(v[i] * result.p + result.q, result.p);
  end;
end;

function abs(f: rational): rational;
begin
  result := frc(abs(f.p), abs(f.q));
end;

function min(f, g: rational): rational;
begin
  if f < g then result := f else result := g;
end;

function max(f, g: rational): rational;
begin
  if f > g then result := f else result := g;
end;

procedure exchange(var a, b: rational);
var
  c: rational;
begin
  c := a;  a := b;  b := c;
end;

procedure order(var a, b: rational); overload;
begin
  if a > b then exchange(a, b);
end;

procedure order(var a, b, c: rational); overload;
begin
  order(a, b);  order(a, c);  order(b, c);
end;

function min(a, b, c: rational): rational;
begin
  order(a, b, c);  result := a;
end;

function mid(a, b, c: rational): rational;
begin
  order(a, b, c);  result := b;
end;

function max(a, b, c: rational): rational;
begin
  order(a, b, c);  result := c;
end;

function censor(l, n, r: rational): rational;
begin
  if n < l then result := l else
  if n > r then result := r else
    result := n;
end;

function fractional(x: numeric; var f: rational; depth: integer; eps: numeric): numeric;
begin
  result := fractional(x, f.p, f.q, depth, eps);
end;

function simplify(f: rational; depth: integer): rational;
begin
  fractional(f.x, result, depth);
end;

function crossprod(f, g: rational): ordinal;
begin
  result := f.p * g.q - f.q * g.p;
end;

procedure fareyseq(n: integer; var f: ratarray);
var
  i: integer;
  o: ordinal;
begin
  if n > 0 then begin
    setlength(f, 2);
    f[0] := 0;  i := 1;  f[i] := frc(1, n);
    while f[i].p <> f[i].q do begin
      i := i + 1;  setlength(f, i + 1);
      o :=(n + f[i-2].q) div f[i-1].q;
      f[i] := frc(o * f[i-1].p - f[i-2].p, o * f[i-1].q - f[i-2].q);
    end;
  end else setlength(f, 0);
end;

function arctan(f: rational): numeric;
begin
  result := atan2(f.p, f.q);
end;


function sqrabs(q: quaternion): numeric;
begin
  // Lagrange: quaternion is prime <=> result is prime
  result := sqr(q.w) + sqr(q.x) + sqr(q.y) + sqr(q.z);
end;

function abs(q: quaternion): numeric;
var
  i: integer;
  a: numeric;
  f: boolean;
begin
  result := abs(q.s);
  i := 1;  f := true;
  while f and (i <= 3) do begin
    a := abs(q.v[i]);
    if a = 0 then i := i + 1 else
      if result = 0 then begin
        result := a;  i := i + 1;
      end else f := false;
  end;
  if (not f) and (i <= 3) then
    result := sqrt(sqrabs(q));
end;

function vectorpart(q: quaternion): quaternion;
begin
  result := q;  result.w := 0;
end;

function scalarpart(q: quaternion): quaternion;
begin
  result := q.w;
end;

function sign(q: quaternion): quaternion;
begin
  if q = 0 then result := q else result := q/abs(q);
end;

function versor(q: quaternion): quaternion;
begin
  q.w := 0;  result := sign(q);
end;


function join(q: quaternion; var r: quaternion): complex;
var
  a: numeric;
begin
  if isnan(q) or isinf(q) then begin
    r := q;
    result := q.f;
  end else begin
    r := vectorpart(q);
    a := abs(r);  if a > 0 then r := r/a;
    result := xiy(q.w, a);
  end;
end;

function split(z: complex; var r: quaternion): quaternion;
begin
  if isnan(z) or isinf(r) then r := z else
  if isnan(r) or isinf(r) then r := z else
    r := z.x + r * z.y;
  result := r;
end;


function sqr(q: quaternion): quaternion;
begin
  result := q * q;
end;

function arg(q: quaternion): quaternion;
var
  a: numeric;
begin
  if q = 0 then result := q else begin
    a := sign(q).s;  q.w := 0;
    if q = 0 then result := q else begin
      a := arccos(a);
      if a = 0
        then result := 0
        else result := a * sign(q);
    end
  end;
end;

function polar(q: quaternion): quaternion;
var
  a: numeric;
  r: complex;
begin
  q.w := 0;
  if q = 0 then result := 1 else begin
    a := abs(q);  r := cis(a);
    result := r.x + r.y * q/a;
  end;
end;

function polar(x: numeric; y: numeric = 0; z: numeric = 0): quaternion;
begin
  result := polar(xyz(x, y, z));
end;

function rotation(q: quaternion): quaternion;
begin
  result := polar(q/2);
end;

function rotation(x: numeric; y: numeric = 0; z: numeric = 0): quaternion;
begin
  result := rotation(xyz(x, y, z));
end;

function angles(q: quaternion): quaternion;
begin
  result := 2 * arg(q);
end;


function exp(q: quaternion): quaternion;
begin
  if q.isscalar then result := exp(q.w) else
  if q.isa then result.a := exp(q.a) else
  if q.isb then result.b := exp(q.b) else
  if q.isc then result.c := exp(q.c) else
    result := exp(q.w) * polar(q);
end;

function ln(q: quaternion): quaternion;
begin
  if q = 0 then result := qtninf else
  if q.isa then result.a := ln(q.a) else
  if q.isb then result.b := ln(q.b) else
  if q.isc then result.c := ln(q.c) else
    result := ln(abs(q)) + arg(q);
end;

function sqrt(q: quaternion): quaternion;
const
  scatter = false;
begin
  if q = 0 then result := 0 else
  if q.isscalar then begin
    result := sqrt(q.f);
    if scatter then
    with result do begin
      x := x/sqrt3;
      y := x;
      z := x;
    end;
  end else
  if q.isa then result.a := sqrt(q.a) else
  if q.isb then result.b := sqrt(q.b) else
  if q.isc then result.c := sqrt(q.c) else
    result := exp(ln(q)/2);
end;


function sin(q: quaternion): quaternion;
begin
  split(sin(join(q, result)), result);
end;

function cos(q: quaternion): quaternion;
begin
  split(cos(join(q, result)), result);
end;

function tan(q: quaternion): quaternion;
begin
  split(tan(join(q, result)), result);
end;

function cot(q: quaternion): quaternion;
begin
  split(cot(join(q, result)), result);
end;


function sinh(q: quaternion): quaternion;
begin
  split(sinh(join(q, result)), result);
end;

function cosh(q: quaternion): quaternion;
begin
  split(cosh(join(q, result)), result);
end;

function tanh(q: quaternion): quaternion;
begin
  split(tanh(join(q, result)), result);
end;

function coth(q: quaternion): quaternion;
begin
  split(coth(join(q, result)), result);
end;


function arcsin(q: quaternion): quaternion;
begin
  split(arcsin(join(q, result)), result);
end;

function arccos(q: quaternion): quaternion;
begin
  split(arccos(join(q, result)), result);
end;

function arctan(q: quaternion): quaternion;
begin
  split(arctan(join(q, result)), result);
end;

function arccot(q: quaternion): quaternion;
begin
  split(arccot(join(q, result)), result);
end;


function arcsinh(q: quaternion): quaternion;
begin
  split(arcsinh(join(q, result)), result);
end;

function arccosh(q: quaternion): quaternion;
begin
  split(arccosh(join(q, result)), result);
end;

function arctanh(q: quaternion): quaternion;
begin
  split(arctanh(join(q, result)), result);
end;

function arccoth(q: quaternion): quaternion;
begin
  split(arccoth(join(q, result)), result);
end;

function Gamma(q: quaternion): quaternion; overload;
begin
  split(Gamma(join(q, result)), result);
end;


function amean(p, q: quaternion; t: numeric = 1/2): quaternion;
begin
  result := p + t * (q - p);
end;

function dotprod(p, q: quaternion): numeric;
begin
  result := p.x * q.x + p.y * q.y + p.z * q.z;
end;

function crossprod(p, q: quaternion): quaternion;
begin
  //result := vectorpart(p) * vectorpart(q);
  result := xyz(
    p.y * q.z - p.z * q.y,
    p.z * q.x - p.x * q.z,
    p.x * q.y - p.y * q.x
  );
end;

function dotcross(p, q: quaternion): quaternion;
begin
  result := dotprod(p, q) + crossprod(p, q);
end;

function intersection(p, q: quaternion): complex;
begin
  result := intersection(p.f, p.g, q.f, q.g);
end;


function qtn(w: numeric; x: numeric = 0; y: numeric = 0; z: numeric = 0): quaternion;
begin
  result.w := w;  result.x := x; result.y := y; result.z := z;
end;

function qtn(u, v: complex): quaternion;
begin
  result := qtn(u.x, u.y, v.x, v.y);
end;

function xyz(x, y, z: numeric): quaternion;
begin
  result := qtn(0, x, y, z);
end;

class operator quaternion.implicit(x: numeric): quaternion;
begin
  result := qtn(x);
end;

class operator quaternion.implicit(z: complex): quaternion;
begin
  result.a := z;
end;

class operator quaternion.positive(q: quaternion): quaternion;
begin
  result := qtn(q.w, q.x, q.y, q.z);
end;

class operator quaternion.negative(q: quaternion): quaternion;
begin
  result := qtn(-q.w, -q.x, -q.y, -q.z);
end;

class operator quaternion.logicalnot(q: quaternion): quaternion;
begin
  result := qtn(q.w, -q.x, -q.y, -q.z);
end;

class operator quaternion.equal(p, q: quaternion): boolean;
begin
  result := (p.w = q.w) and (p.x = q.x) and (p.y = q.y) and (p.z = q.z);
end;

class operator quaternion.notequal(p, q: quaternion): boolean;
begin
  result := not (p = q);
end;

class operator quaternion.add(p, q: quaternion): quaternion;
begin
  result := qtn(p.w + q.w, p.x + q.x, p.y + q.y, p.z + q.z);
end;

class operator quaternion.add(q: quaternion; x: numeric): quaternion;
begin
  result := q;  result.w := result.w + x;
end;

class operator quaternion.add(x: numeric; q: quaternion): quaternion;
begin
  result := q + x;
end;

class operator quaternion.subtract(p, q: quaternion): quaternion;
begin
  result := -q + p;
end;

class operator quaternion.subtract(q: quaternion; x: numeric): quaternion;
begin
  result := -x + q;
end;

class operator quaternion.subtract(x: numeric; q: quaternion): quaternion;
begin
  result := -q + x;
end;

class operator quaternion.multiply(x: numeric; q: quaternion): quaternion;
begin
  result := qtn(x * q.w, x * q.x, x * q.y, x * q.z);
end;

class operator quaternion.multiply(q: quaternion; x: numeric): quaternion;
begin
  result := x * q;
end;

class operator quaternion.multiply(p, q: quaternion): quaternion;
begin
  result := qtn(
    p.w*q.w - p.x*q.x - p.y*q.y - p.z*q.z,
    p.w*q.x + p.x*q.w + p.y*q.z - p.z*q.y,
    p.w*q.y - p.x*q.z + p.y*q.w + p.z*q.x,
    p.w*q.z + p.x*q.y - p.y*q.x + p.z*q.w
  );
end;

class operator quaternion.divide(q: quaternion; x: numeric): quaternion;
begin
  if x = 0
    then result := qtninf
    else result := qtn(q.w/x, q.x/x, q.y/x, q.z/x);
end;

class operator quaternion.divide(p, q: quaternion): quaternion;
begin
  result := p * (not q) / sqrabs(q);
end;

class operator quaternion.divide(x: numeric; q: quaternion): quaternion;
begin
  result := x * (not q) / sqrabs(q);
end;

class operator quaternion.bitwiseand(p, q: quaternion): quaternion;
begin
  result := (q * p) * (not q);
end;

class operator quaternion.leftshift(p, q: quaternion): quaternion;
var
  a: numeric;
begin
  a := abs(q);  if a = 0 then q := 1 else q := q/a;
  result := p and q;  result.w := p.w;
end;

class operator quaternion.rightshift(p, q: quaternion): quaternion;
begin
  result := p shl (not q);
end;


function quaternion.getisscalar: boolean;
begin
  result := vectorpart(self) = 0;
end;

procedure quaternion.setisscalar(f: boolean);
begin
  if f then begin
    x := 0; y := 0; z := 0;
  end;
end;

function quaternion.getisvector: boolean;
begin
  result := w = 0;
end;

procedure quaternion.setisvector(f: boolean);
begin
  if f then w := 0;
end;

function quaternion.geta: complex;
begin
  result := xiy(w, x);
end;

procedure quaternion.seta(z: complex);
begin
  self := qtn(z.x, z.y, 0, 0);
end;

function quaternion.getb: complex;
begin
  result := xiy(w, y);
end;

procedure quaternion.setb(z: complex);
begin
  self := qtn(z.x, 0, z.y, 0);
end;

function quaternion.getc: complex;
begin
  result := xiy(w, z);
end;

procedure quaternion.setc(z: complex);
begin
  self := qtn(z.x, 0, 0, z.y);
end;

function quaternion.getisa: boolean;
begin
  result := (y = 0) and (z = 0);
end;

function quaternion.getisb: boolean;
begin
  result := (x = 0) and (z = 0);
end;

function quaternion.getisc: boolean;
begin
  result := (x = 0) and (y = 0);
end;

procedure exchange(var a, b: quaternion);
var
  c: quaternion;
begin
  c := a;  a := b;  b := c;
end;

function power(q: quaternion; x: numeric): quaternion;
begin
  if (q = 1) or  (x = 0) then result := 1 else
  if (q = 0) and (x > 0) then result := q else
  if x =   1  then result := q else
  if x =  -1  then result := 1/q else
  if x =   2  then result := sqr(q) else
  if x =  -2  then result := 1/sqr(q) else
  if x =  1/2 then result := sqrt(q) else
  if x = -1/2 then result := 1/sqrt(q) else
  if q.isa  then result.a := power(q.a, x) else
  if q.isb  then result.b := power(q.b, x) else
  if q.isc  then result.c := power(q.c, x) else
    result := exp(ln(q) * x);
end;

function power(q: quaternion; z: complex): quaternion;
begin
  if z.isreal then result := power(q, z.x) else
    split(power(join(q, result), z), result);
end;

function power(p, q: quaternion): quaternion;
begin
  if q.isscalar then result := power(p, q.w) else
  if q.isa then result := power(p, q.a) else
    result := exp(ln(p) * q);
end;


function mandelbrot(c: quaternion; x: numeric = 2; m: integer = 1000): integer;
var
  q, p: quaternion;
begin
  c.w := 0;  result := m;   q := c;
  while (result > 0) and (sqrabs(q) <= 4) do begin
    p := q;  q := power(q, x) + c;
    if q = p
      then result := -result
      else result := result - 1;
  end;
end;


function yprtoqtn(q: quaternion): quaternion;
var
  y, p, r: complex;
begin
  y := cis(q.x/2);  p := cis(q.y/2);  r := cis(q.z/2);
  result := qtn(
    y.x * p.x * r.x + y.y * p.y * r.y,
    y.x * p.x * r.y - y.y * p.y * r.x,
    y.y * p.x * r.y + y.x * p.y * r.x,
    y.y * p.x * r.x - y.x * p.y * r.y
  );
end;

function yprtoqtn(yaw, pitch, roll: numeric): quaternion;
begin
  result := yprtoqtn(xyz(yaw, pitch, roll) * rad);
end;

function qtntoypr(q: quaternion): quaternion;
const
  gimbal_lock: numeric = 1e-19;
begin
  if q = 0 then result := 0 else begin
    result.y := 2*(q.w * q.y - q.z * q.x);
    if (1 - abs(result.y)) < gimbal_lock
      then result := xyz(
        0,
        psign(result.y, π/2),
        2*arg(q.w, q.x)
      )
      else result := xyz(
        arg(1 - 2*(sqr(q.y) + sqr(q.z)), 2*(q.w * q.z + q.x * q.y)),
        arcsin(result.y),
        arg(1 - 2*(sqr(q.x) + sqr(q.y)), 2*(q.w * q.x + q.y * q.z))
      )
  end;
end;

function qtntoypr(q: quaternion; var yaw, pitch, roll: numeric): quaternion;
begin
  result := qtntoypr(q);
  yaw := deg(result.x);  pitch := deg(result.y);  roll := deg(result.z);
end;


var seed: natural = 0;

procedure entropy(var seed: natural);
function rdtsc: natural; assembler; // read time-stamp counter
asm dw 310fh end; // no mnemonic, use opcode
begin
  seed := rdtsc;  seed := mmixrng(seed);
end;

function mmixrng(var seed: natural): natural; // mixed Knuth MMIX lcprng & Marsaglia xorshift
const
  a = $5851F42D4C957F2D; // b = 1/a = multinv(a)
  c = $14057B7EF767814F; // d = -bc = not (b*c) + 1
  x = 13; y = 9; z = 15;
begin
  seed := seed * a + c; result := seed;
  result := result xor (result shl x);
  result := result xor (result shr y);
  result := result xor (result shl z);
end;

function mmixrng: natural;
begin
  result := mmixrng(seed);
end;

function ordrandom(var seed: natural; n: ordinal): ordinal;
var
  f: natural;
begin
  if n > 1 then begin
    // himul(f, n, mmixrng(seed));  result := f;
    f := basediv(n); // f = 2⁶⁴ div n
    repeat
      result := mmixrng(seed) div f;
    until result < n;
  end else result := 0;
end;

function ordrandom(n: ordinal): ordinal; overload;
begin
  result := ordrandom(seed, n);
end;

function ordrange(var seed: natural; n1, n2: ordinal): ordinal;
begin
  order(n1, n2);
  result := n1 + ordrandom(seed, n2 - n1 + 1);
end;

function ordrange(n1, n2: ordinal): ordinal;
begin
  result := ordrange(seed, n1, n2);
end;

function ratrandom(var seed: natural; range: ordinal): rational; overload;
begin
  result := frc(ordrandom(seed, range), range);
end;

function ratrandom(range: ordinal = 510510): rational; overload;
begin
  result := ratrandom(seed, range);
end;

function unidev(var seed: natural; x: numeric = 1): numeric;
begin
  repeat
    result := mmixrng(seed) * nateps;
  until result < 1;
  result := result * x;
end;

function unidev(x: numeric = 1): numeric;
begin
  result := unidev(seed, x);
end;

function unidev(var seed: natural; a, b: numeric): numeric;
begin
  result := a + unidev(seed, b - a);
end;

function unidev(a, b: numeric): numeric;
begin
  result := unidev(seed, a, b);
end;

function expdev(var seed: natural; lambda: numeric = 1): numeric;
begin
  if lambda > 0
    then result := lnone(unidev(seed)) / lambda
    else result := 0;
end;

function expdev(lambda: numeric = 1): numeric;
begin
  result := expdev(seed, lambda);
end;

function geomdev(var seed: natural; p: numeric): ordinal;
begin
  if (0 < p) and (p < 1)
    then result := trunc(expdev(seed, lnone(p)))
    else result := 0;
end;

function geomdev(p: numeric): ordinal;
begin
  result := geomdev(seed, p);
end;

function polardev(var seed: natural; r: numeric): complex;
begin
  result := sqrt(abs(r)) * cis(unidev(seed, τ));
end;

function polardev(r: numeric = 1): complex;
begin
  result := polardev(seed, r);
end;

function gaussdev(var seed: natural): complex;
begin
  result := polardev(seed, expdev(seed, 1/2));
end;

function gaussdev: complex;
begin
  result := gaussdev(seed);
end;

function normaldev(var seed: natural; mu, sigma: complex; angle: numeric): complex;
begin
  if sigma = 0 then result := mu else begin
    result := gaussdev(seed);
    result := xiy(result.x * sigma.x, result.y * sigma.y) * cis(angle) + mu;
  end;
end;

function normaldev(mu, sigma: complex; angle: numeric): complex;
begin
  result := normaldev(seed, mu, sigma, angle);
end;

function erlangdev(var seed: natural; k: integer; lambda: numeric = 1): numeric;
begin
  if lambda > 0 then begin
    result := 1;
    while k > 0 do begin
      result := result - unidev(seed, result);
      k := k - 1;
    end;
    result := -ln(result) / lambda;
  end else result := 0;
end;

function erlangdev(k: integer; lambda: numeric = 1): numeric;
begin
  result := erlangdev(seed, k, lambda);
end;

function chi2dev(var seed: natural; nu: numeric): numeric;
var
  n: integer;
begin
  if nu > 0 then begin
    if frac(nu) = 0 then begin
      n := trunc(nu);
      result := 2 * erlangdev(seed, n div 2);
      if odd(n) then result := result + sqr(gaussdev(seed).x);
    end else result := 2 * gammadev(seed, nu/2);
  end else result := 0;
end;

function chi2dev(nu: numeric): numeric;
begin
  result := chi2dev(seed, nu);
end;

function studenttdev(var seed: natural; nu: numeric): numeric;
begin
  if nu > 0 then begin
    repeat result := chi2dev(seed, nu); until result > 0;
    result := gaussdev(seed).x * sqrt(nu/result);
  end else result := 0;
end;

function studenttdev(nu: numeric): numeric;
begin
  result := studenttdev(seed, nu);
end;

function gammadev(var seed: natural; a: numeric): numeric;
var
  x: numeric;
  f: boolean;
begin
  if a > 0 then begin
    result := erlangdev(seed, trunc(a));
    a := frac(a);
    if a > 0 then begin
      if a = 1/2 then x := chi2dev(seed, 1)/2 else
      repeat
        if unidev(seed, ℮ + a) < ℮ then begin
          x := power(unidev(seed), 1/a);
          f := toss(seed, exp(-x));
        end else begin
          x := 1 + expdev(seed);
          f := toss(seed, power(x, a - 1));
        end;
      until f;
      result := result + x;
    end;
  end else result := 0;
end;

function gammadev(a: numeric): numeric;
begin
  result := gammadev(seed, a);
end;

function poissondev(var seed: natural; lambda: numeric): ordinal;
var
  p, l: numeric;
begin
  if lambda > 0 then begin
    l := exp(-lambda);  result := -1;  p := 1;
    repeat
      result := result + 1;  p := unidev(seed, p);
    until p <= l;
  end else result := 0;
end;

function poissondev(lambda: numeric): ordinal;
begin
  result := poissondev(seed, lambda);
end;

function toss(var seed: natural; p: numeric = 1/2): boolean;
begin
  result := unidev(seed) < p;
end;

function toss(p: numeric = 1/2): boolean;
begin
  result := toss(seed, p);
end;

function berndev(var seed: natural; n, k: ordinal): boolean;
begin
  result := (n >= 0) and (k >= 0);
  if result and (n > k) then result := ordrandom(seed, n) < k;
end;

function berndev(n, k: ordinal): boolean;
begin
  result := berndev(seed, n, k)
end;

function bindev(var seed: natural; n: ordinal; p: numeric = 1/2): ordinal;
const
  limit = 64;
var
  i: integer;
  mean, deviate: numeric;
begin
  result := 0;
  if (n > 0) and (p > 0) then
  if p >= 1 then result := n else
  begin
    if n > limit then begin
      mean := n * p;  deviate := sqrt(mean * (1 - p));
      repeat
        result := round(normaldev(seed, mean, deviate).x);
      until (0 <= result) and (result <= n);
    end else begin
      for i := 1 to n do if toss(seed, p) then result := result + 1;
    end;
  end;
end;

function bindev(n: ordinal; p: numeric = 1/2): ordinal;
begin
  result := bindev(seed, n, p);
end;

function benforddev(var seed: natural; m, n: ordinal): ordinal;
begin
  order(m, n);
  if m > 0
    then result := censor(m, trunc(exp10(unidev(seed, log10(m), log10(n + 1)))), n)
    else result := 0;
end;

function benforddev(m, n: ordinal): ordinal;
begin
  result := benforddev(seed, m, n);
end;

procedure ordstir(var seed: natural; var x: array of ordinal; n: integer);
var
  i: integer;
begin
  if n < 0
    then n := high(x) - 1
    else n := min(n, length(x)) - low(x) - 1;
  for i := low(x) to low(x) + n do
    exchange(x[i], x[ordrange(seed, i, high(x))]);
end;

procedure ordstir(var x: array of ordinal; n: integer);
begin
  ordstir(seed, x, n);
end;

procedure ordfill(var seed: natural; var x: array of ordinal; ordered: boolean = false);
var
  i, j, k: integer;
begin
  k := 1;
  for i := low(x) to high(x) do begin
    if ordered
      then j := i
      else j := ordrange(seed, low(x), i);
    x[i] := x[j];  x[j] := k;
    k := k + 1;
  end;
end;

procedure ordfill(var x: array of ordinal; ordered: boolean = false);
begin
  ordfill(seed, x, ordered);
end;

procedure ordsample(var seed: natural; var x: array of ordinal; n: integer);
var
  i, j, k: integer;
function iter(var l: integer): integer; // l++
begin
  result := l;  l := l + 1;
end;
begin
  i := 0;  j := low(x);  k := length(x) - j;  if n < k then n := k;
  while j < k do if berndev(seed, n - iter(i), k - j) then x[iter(j)] := i;
end;

procedure ordsample(var x: array of ordinal; n: integer);
begin
  ordsample(seed, x, n);
end;

procedure ordsort(var x: array of ordinal; lpart, rpart: integer);
const
  treshold = 16;
var
  l, r: integer;
  p: ordinal;
begin
  if (rpart - lpart) > treshold then begin
    // Quick Sort (partition-exchange),
    // Sir Charles Antony Richard Hoare (1960)
    l := lpart;  r := rpart;
    p := x[ordrange(l, r)]; // random pivot: 1.386 * O(n ln n)
    repeat
      while x[l] < p do l := l + 1;
      while x[r] > p do r := r - 1;
      if l <= r then begin
        exchange(x[l], x[r]);
        l := l + 1;  r := r - 1;
      end;
    until l > r;
    if r > lpart then ordsort(x, lpart, r);
    if l < rpart then ordsort(x, l, rpart);
  end else begin
    // Insertion Sort is better for short arrays
    for r := lpart + 1 to rpart do begin
      l := r;
      p := x[l];
      while (l > lpart) and (x[l - 1] > p) do begin
        x[l] := x[l - 1];  l := l - 1;
      end;
      x[l] := p;
    end;
  end;
end;

procedure ordsort(var x: array of ordinal);
begin
  if length(x) > 1 then ordsort(x, low(x), high(x));
end;


procedure numsort(var x: array of numeric; lpart, rpart: integer);
const
  treshold = 16;
var
  l, r: integer;
  p: numeric;
begin
  if (rpart - lpart) > treshold then begin
    l := lpart;  r := rpart;
    p := x[ordrange(l, r)];
    repeat
      while x[l] < p do l := l + 1;
      while x[r] > p do r := r - 1;
      if l <= r then begin
        exchange(x[l], x[r]);
        l := l + 1;  r := r - 1;
      end;
    until l > r;
    if r > lpart then numsort(x, lpart, r);
    if l < rpart then numsort(x, l, rpart);
  end else begin
    for r := lpart + 1 to rpart do begin
      l := r;
      p := x[l];
      while (l > lpart) and (x[l - 1] > p) do begin
        x[l] := x[l - 1];  l := l - 1;
      end;
      x[l] := p;
    end;
  end;
end;

procedure numsort(var x: array of numeric);
begin
  if length(x) > 1 then numsort(x, low(x), high(x));
end;



function zxf(n: natural): zxfloat;

begin

  result.b := n and 255;  n := n shr 8;
  result.c := n and 255;  n := n shr 8;

  result.d := n and 255;  n := n shr 8;

  result.e := n and 255;  n := n shr 8;

  result.a := n and 255;  n := n shr 8;

end;


function zxf(a, e, d, c, b: byte): zxfloat;

begin

  result.a := a;

  result.b := b;
  result.c := c;
  result.d := d;
  result.e := e;
end;


function zxfloat.getm: numeric;

begin

  if a = 0 then begin
    result := c; result := $100*result + d;
    if e = $ff then result := result - $10000;
  end else begin

    result := (((b/$100 + c)/$100 + d)/$100 + (e or $80))/$100;

    if e and $80 <> 0 then result := -result;
  end;
end;


function zxfloat.gete: integer;

begin

  result := a;

  if result > 0 then result := result - $80;
end;


function zxfloat.getx: numeric;

begin

  result := mantissa * power(2, exponent);

end;


procedure zxfloat.setx(r: numeric);

var

  l: integer;

  n: natural;

  s: boolean;

begin

  if isnan(r) or (r = 0) then self := zxzero else

  if r > zxmax.x then self := zxmax else

  if r < zxmin.x then self := zxmin else

  begin

    s := r < 0;  r := abs(r);  l := ilog2(r) + $80;

    while r >= 1 do begin

      r := r/2;  l := l + 1;
    end;

    n := round(r * $100000000);

    if n >= $100000000 then begin

      n := n div 2;  l := l + 1;

    end;

    if l > 0 then begin

      a := l;

      b := n and $ff;  n := n shr 8;

      c := n and $ff;  n := n shr 8;

      d := n and $ff;  n := n shr 8;

      e := n and $7f;  n := n shr 7;

      if s then e := e or $80;

      r := x;

      if (frac(r) = 0) and (-$10000 <= r) and (r < $10000) then begin

        l := trunc(r);
        a := 0;  b := 0;
        d := l and $ff;  l := l shr 8;
        c := l and $ff;  l := l shr 8;

        e := l and $ff;  l := l shr 8;

      end;

    end else self := zxzero;

  end;

end;


class operator zxfloat.implicit(z: zxfloat): numeric;

begin

  result := z.x;

end;


class operator zxfloat.implicit(x: numeric): zxfloat;
begin
  result.x := x;
end;


class operator zxfloat.equal(f, g: zxfloat): boolean;
begin
  result := f.x = g.x;
end;


class operator zxfloat.notequal(f, g: zxfloat): boolean;
begin
  result := f.x <> g.x;
end;


class operator zxfloat.lessthan(f, g: zxfloat): boolean;
begin
  result := f.x < g.x;
end;


class operator zxfloat.greaterthan(f, g: zxfloat): boolean;
begin
  result := f.x > g.x;
end;


class operator zxfloat.lessthanorequal(f, g: zxfloat): boolean;
begin
  result := f.x <= g.x;
end;


class operator zxfloat.greaterthanorequal(f, g: zxfloat): boolean;
begin
  result := f.x >= g.x;
end;


class operator zxfloat.positive(f: zxfloat): zxfloat;
begin
  result := +f.x;
end;


class operator zxfloat.negative(f: zxfloat): zxfloat;
begin
  result := -f.x;
end;

class operator zxfloat.add(f, g: zxfloat): zxfloat;
begin
  result := f.x + g.x;
end;

class operator zxfloat.subtract(f, g: zxfloat): zxfloat;
begin
  result := f.x - g.x;
end;

class operator zxfloat.multiply(f, g: zxfloat): zxfloat;
begin
  result := f.x * g.x;
end;

class operator zxfloat.divide(f, g: zxfloat): zxfloat;
begin
  result := f.x / g.x;
end;


class operator zxcomplex.implicit(x: numeric): zxcomplex;

begin

  result.x := x;  result.y := 0;
end;


class operator zxcomplex.implicit(c: zxcomplex): complex;
begin
  result.x := c.x;  result.y := c.y;
end;


class operator zxcomplex.implicit(z: complex): zxcomplex;
begin
  result.x := z.x;  result.y := z.y;
end;




procedure init;

{$IFDEF WIRTH}
var
  a: numeric;
{$ELSE}
function ldln2: numeric;
asm
  FLDLN2; FWAIT;
end;
function ldlg10: numeric;
asm
  FLDL2T; FWAIT;
end;
{$ENDIF}
begin
  nateps := power(2, -8 * sizeof(natural));
  entropy(seed);
  calctinies;
  {$IFDEF WIRTH}
  a := 1; repeat sqrt2 := a; a := a/2 + 1/a; until sqrt2 = a;
  {$ELSE}
  sqrt2 := sqrt(2);
  {$ENDIF}
  sqrth := sqrt2/2;
  sqrt3 := sqrt(3);
  sqrt5 := sqrt(5);
  π := pi;  τ := 2*π;
  sqrtπ := sqrt(π);
  φ := amean(sqrt5, 1);
  {$IFDEF WIRTH}
  if ln2 = 0 then ln2 := lniter(1, 1/9) / 3; // calculate 6 arccoth(3) = 3 ln(2)
  if lg10 = 0 then lg10 := lniter(1/9) / ln2 + 3;
  if ℮ = 0 then taylorexp(1, a, a, a, a, ℮);
  {$ELSE}
  ln2 := ldln2;
  lg10 := ldlg10;
  ℮ := exp2(1/ln2);
  {$ENDIF}
end;


procedure done;
begin
end;


initialization


init;


finalization


done;


end.

