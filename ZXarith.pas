unit ZXarith;


interface


uses Blaise;


type
  zxfloat = packed record
    a, e, d, c, b: byte;
  private
    function getx: numeric;
    procedure setx(x: numeric);
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
    property value: numeric read getx write setx;
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
    class operator logicalnot(c: zxcomplex): zxcomplex;

    class operator equal(u, v: zxcomplex): boolean;

    class operator notequal(u, v: zxcomplex): boolean;
end;


implementation


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

procedure zxfloat.setx(x: numeric);
var
  l: integer;
  n: natural;
  s: boolean;
begin
  if isnan(x) or (x = 0) then self := zxzero else
  if x > zxmax.value then self := zxmax else
  if x < zxmin.value then self := zxmin else
  begin
    s := x < 0;  x := abs(x);  l := ilog2(x) + $80;
    while x >= 1 do begin
      x := x/2;  l := l + 1;
    end;
    n := round(x * $100000000);
    if n >= $100000000 then begin
      n := n div 2;  l := l + 1;
    end;
    if l > 0 then begin
      a := l;
      b := n and $ff;  n := n shr 8;
      c := n and $ff;  n := n shr 8;
      d := n and $ff;  n := n shr 8;
      e := n and $7f;
      if s then e := e or $80;
      x := value;
      if (frac(x) = 0) and (-$10000 <= x) and (x < $10000) then begin
        l := trunc(x);
        a := 0;  b := 0;
        d := l and $ff;  l := l shr 8;
        c := l and $ff;  l := l shr 8;
        e := l and $ff;
      end;
    end else self := zxzero;
  end;
end;

class operator zxfloat.implicit(z: zxfloat): numeric;
begin
  result := z.value;
end;

class operator zxfloat.implicit(x: numeric): zxfloat;
begin
  result.value := x;
end;

class operator zxfloat.equal(f, g: zxfloat): boolean;
begin
  result := f.value = g.value;
end;


class operator zxfloat.notequal(f, g: zxfloat): boolean;
begin
  result := f.value <> g.value;
end;


class operator zxfloat.lessthan(f, g: zxfloat): boolean;
begin
  result := f.value < g.value;
end;


class operator zxfloat.greaterthan(f, g: zxfloat): boolean;
begin
  result := f.value > g.value;
end;

class operator zxfloat.lessthanorequal(f, g: zxfloat): boolean;
begin
  result := f.value <= g.value;
end;

class operator zxfloat.greaterthanorequal(f, g: zxfloat): boolean;
begin
  result := f.value >= g.value;
end;

class operator zxfloat.positive(f: zxfloat): zxfloat;
begin
  result := +f.value;
end;

class operator zxfloat.negative(f: zxfloat): zxfloat;
begin
  result := -f.value;
end;

class operator zxfloat.add(f, g: zxfloat): zxfloat;
begin
  result := f.value + g.value;
end;

class operator zxfloat.subtract(f, g: zxfloat): zxfloat;
begin
  result := f.value - g.value;
end;

class operator zxfloat.multiply(f, g: zxfloat): zxfloat;
begin
  result := f.value * g.value;
end;

class operator zxfloat.divide(f, g: zxfloat): zxfloat;
begin
  result := f.value / g.value;
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


class operator zxcomplex.logicalnot(c: zxcomplex): zxcomplex;

begin
  result.x := +c.x;  result.y := -c.y;
end;


class operator zxcomplex.equal(u, v: zxcomplex): boolean;

begin

  result := (u.x = v.x) and (u.y = v.y);
end;


class operator zxcomplex.notequal(u, v: zxcomplex): boolean;
begin
  result := not (u = v);
end;



function zxf(n: natural): zxfloat;

begin
  result.b := n and 255;  n := n shr 8;
  result.c := n and 255;  n := n shr 8;
  result.d := n and 255;  n := n shr 8;
  result.e := n and 255;  n := n shr 8;
  result.a := n and 255;
end;

function zxf(a, e, d, c, b: byte): zxfloat;
begin
  result.a := a;
  result.b := b;
  result.c := c;
  result.d := d;
  result.e := e;
end;



end.
