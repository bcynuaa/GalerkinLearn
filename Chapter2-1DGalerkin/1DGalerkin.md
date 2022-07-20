<h1>一维问题的有限元伽辽金解法</h1>

[toc]

# 0.问题的引入

虽然说Galerkin法非常好用，但在面对更加复杂的问题时未免捉襟见肘，精度不够。在整个计算域上用一个函数表达式来进行描述固然很方便，但带来的问题会有精度上的不足，而且解得最终效果也很依赖于基函数的选取。

那一个很自然的想法：**能否将整个区域分成小块，在小块上进行一些基函数近似。**

为了讨论这个问题，我们设想一个线性微分算子$A$（如$\frac{d^2}{dx^2}+3\frac{d}{dx}$），对于方程$Au=f$而言，我们展开如下的讨论。

# 1.直觉中的FDM解法

早在学习高等数学的时候，我就想着人们完全可以用差分来替代微分以获取“精度不那么高的解”。事实是在计算方法这门课程中，也是这么做的。大概的方法就是我们可以根据极限的定义，用差商代替微商：

$$
\lim_{h\to 0}\frac{y(x+h)-y(x)}{h}=y'(x)
$$

实施过程中我们不可能把$h$真的控制为$0$，因此一般只选取较小的步长。

我们来考虑一个横跨$[a,b]$区间的微分方程$y'=f$的数值解：

首先将区间划分为$n$份，步长$h=\frac{b-a}{n}$，接下来将划分好的区间的$n+a$个节点上标上编号，依次记作$x_0,x_1,x_2\dots,x_n$，有$x_i=x_0+ih$，并将其对应的函数值$y$记作$y_i$。

接下来，我们用泰勒展开，将$y_{i+1}$用$y_i$表示：

$$
y_{i+1}=y_i + y'(x_i)h + \frac{y''(x_i)}{2}h^2 +\frac{y'''(x_i)}{6}h^3\dots
$$

因为在划分时$h$是一个小量，因此高阶项可以近似忽略不计，我们将上式近似表示为：

$$
\frac{y_{i+1}-y_{i}}{h}\approx y'(x_i)=f(x_i)
$$

因此在$y_i$已知的情况下，可以递推得到：

$$
y_{i+1}=y_{i}+f(x_i)h
$$

我们可以来思考一个简单的例子：

$$
y'=x^3 e^x\quad y(0)=0\quad y(3)=?
$$

我们用前述的方法（因为这是很直给的一步一步向前迭代计算的方法，因此称为前进欧拉法）进行测试，如代码`Chapter2-1DGalerkin//code//Euler.jl`所述：

```julia
# Chapter2-1DGalerkin//code//Euler.jl
function EulerForward(a, b, y0, f, devide=1000)
    h = (b - a) / devide;
    x = collect(a: h: b);
    y = zeros(devide + 1);
    y[1] = y0;
    for i = 2 : devide+1
        y[i] = y[i-1] + f(x[i]) * h;
    end
    return x, y;
end
```

显然这个方程是有显示解的，我们表示出来并和上述数值解进行比照：

$$
y(x) = (x^3-3x^2+6x-6)e^{x}+6
$$

将比较程序列出，如下：

```julia
# Chapter2-1DGalerkin//code//test.ipynb
using Plots;
include("Euler.jl");

a = 0.;
b = 3.;
y0 = 0;
devide = 300;
h = (b-a)/devide

function f_origin(x)
    y = (x^3 - 3*x^2 + 6*x - 6) * exp(x) + 6;
    return y;
end

function f(x)
    y = exp(x) * x^3;
    return y;
end

xreal = collect(a: h: b);
yreal = f_origin.(xreal);

xnumer, ynumer = EulerForward(a, b, y0, f, devide)

plot(xreal, yreal, lw=10, label="real solution", alpha=0.5)
plot!(xnumer, ynumer, lw=3, label="numerical solution")
savefig("..//image//fdm_example.png")
```

![真实解与解析解的对比 image//fdm_example.png](image//fdm_example.png)

可以发现两个解是相似的。虽然到后面因为$e^x$形式函数增长很快，两条曲线逐渐有差别，但在有限范围内，两条函数曲线吻合得比较好。

# 2.我们能用Galerkin法求解上述问题吗？

FDM已经提供了一个相当直观而且简单有效的微分方程解法了，但我们的重点还是在Galerkin上，我们关注如果用伽辽金法的思想，应该怎么处理这个问题。

如果是一般的伽辽金法，我们该如何操作呢？显然这个时候如果把眼光放在全局，在$[0,3]$区间上进行全局的伽辽金法求解，那么解可能是惨不忍睹的。不妨来试一下：

假设一系列的基函数，取满足条件边界条件的三次函数形式：

$$
\tilde{y} = a_1 x + a_2 x^2 + a^3 x^3
$$

按照前文所述的，记：

$$
H(F) = \int_0^3 F(x)dx
$$

记$f(x)=x^3 e^x$，$\phi_j=x^j,j=1,2,3$。则可以列出方程：

$$
\begin{bmatrix}
H(\phi_1 \phi_1') & H(\phi_1 \phi_2') & H(\phi_1 \phi_3')\\
H(\phi_2 \phi_1') & H(\phi_2 \phi_2') & H(\phi_2 \phi_3')\\
H(\phi_3 \phi_1') & H(\phi_3 \phi_2') & H(\phi_3 \phi_3')
\end{bmatrix}
\begin{bmatrix}
a_1\\a_2\\a_3
\end{bmatrix}=
\begin{bmatrix}
H(f\phi_1)\\H(f\phi_2)\\H(f\phi_3)
\end{bmatrix}
$$

计算过程在`1DGalerkin.nb`中，计算结果如下图所示：

![为什么不用Galerkin image//NoGalerkin.png](image//NoGalerkin.png)

显然Galerkin解和真实解相去甚远——因此我们用一个补救的措施：将区间划分成小段，在每个小段上使用Galerkin法。

# 3.小段上的伽辽金法研究

让我们来考察一个微小的区间$[x_k,x_{k+1}]$，记步长为$h_k=x_{k+1}-x_k$。假定在$x_k$点上，方程$Au=f$的解为$u_k$，在$x_{k+1}上解为$$u_{k+1}$，这样的一个区间被称作一个一维单元。