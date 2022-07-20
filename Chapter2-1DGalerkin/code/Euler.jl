"""
# 前进欧拉解法

EulerForward(a, b, y0, f, devide=1000::Int) -> Array

用于解决形如下述的方程

```math
y'(x)=f(x)
```

其中在区间`[a,b]`上，分成devide个部分，初值为y0的有限差分前进欧拉解法。
"""
function EulerForward(a, b, y0, f, devide=1000::Int)
    h = (b - a) / devide;
    x = collect(a: h: b);
    y = zeros(devide + 1);
    y[1] = y0;
    for i = 2 : devide+1
        y[i] = y[i-1] + f(x[i]) * h;
    end
    return x, y;
end