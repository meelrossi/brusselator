function [c,s] = givensrotation(a,b)
  if b != 0
    if abs(b) > abs(a)
      r = a / b;
      s = 1 / sqrt(1 + r^2);
      c = s*r;
    else
      r = b / a;
      c = 1 / sqrt(1 + r^2);
      s = c*r;
    end
  else
    c = 1;
    s = 0;
  end

end