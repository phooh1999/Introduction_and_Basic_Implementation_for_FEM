# 0 引言

> 希望能吃透有限元算法的构造，尤其是通用结构，而不是case by case处理问题。以及尽可能自己把有限元程序写出来，整个程序自己写出来和直接看参考程序是完全不同的体验。

这里简单记录课程的前一半学习内容，包括
1. 一维二阶椭圆方程
2. 二维二阶椭圆方程
3. 二维二阶抛物方程
4. 二维线弹性方程（狄利克雷边界条件）

# Chapter 1: Finite Elements for 1D second order elliptic equation

## 1.1 Weak/Galerkin formulation

一维稳态线性二阶椭圆，狄利克雷条件

梳理写课程作业程序前的思路，主要目的在于写出程序，因此涉及到的更深入的数学知识如Sobolev Spaces等暂且搁置。

Solve
$$
\begin{gather}
-\frac{\mathrm{d}}{\mathrm{d}x}\left( c\left( x \right) \frac{\mathrm{d}u\left( x \right)}{\mathrm{d}x} \right) =f\left( x \right) ,\ a<x<b
\nonumber
\\
u\left( a \right) =g_a, \  u\left( b \right) =g_b
\nonumber
\end{gather} 
$$
for $u(x)$

两边同时乘以 $v(x)$ 并积分（**注意这里的 $v(x)$ 是可以任意选择的，而 $u(x)$ 由实际的物理决定**）
$$
-\int_a^b{\frac{\mathrm{d}}{\mathrm{d}x}\left( c\left( x \right) \frac{\mathrm{d}u\left( x \right)}{\mathrm{d}x} \right) v\left( x \right) \mathrm{d}x=\int_a^b{f\left( x \right) v\left( x \right) \mathrm{d}x}}
$$
再分部积分
$$
-c\left( b \right) u^\prime\left( b \right) v\left( b \right) +c\left( a \right) u^\prime\left( a \right) v\left( a \right) +\int_a^b{cu^\prime v^\prime\mathrm{d}x=\int_a^b{fv\mathrm{d}x}} 
$$
选取 $v(x)$ 使得 $v(a)=v(b)=0$，因此有
$$
\int_a^b{cu^\prime v^\prime\mathrm{d}x=\int_a^b{fv\mathrm{d}x}}
$$
## 1.2 FE Space/discretization

在 $[a,b]$ 上划分均匀网格 $h=(b-a)/N$
$$
\begin{gather}
x_i=a+\left( i-1 \right) h\,\, \left( i=1,\cdots ,N+1 \right) 
\nonumber
\\
E_n=\left[ x_n,x_{n+1} \right] \,\,\left( n=1,\cdots ,N \right) 
\nonumber
\end{gather}
$$
使用分片线性函数
$$
\phi _j\left( x_i \right) =\delta _{ij}\,\quad i,j=1,\cdots ,N+1
$$
在有限元空间下
$$
\begin{gather}
\int_a^b{cu^\prime_h v^\prime_h\mathrm{d}x=\int_a^b{fv_h\mathrm{d}x}}
\nonumber
\\
u_h=\sum_{j=1}^{N+1}{u_j\phi _j}
\nonumber
\end{gather}
$$
将test function选为同样的分片线性函数 $v_h=\phi_i \  \left( i=1,\cdots,N+1 \right)$ 
> 因为 $u_j \ (j=1,\cdots,N+1)$ 共 $(N+1)$ 个未知数需要求解，所以选择 $(N+1)$ 个test function测试，最后形成方程组联立求解


$$
\begin{aligned}
&\int_a^b{c\left( \sum_{j=1}^{N+1}{u_j\phi _j} \right) ^{\prime}\phi _{i}^{\prime}\mathrm{d}x=\int_a^b{f\phi _{i}\mathrm{d}x}}
,\quad i=1,\cdots,N+1
\\
\Rightarrow &\sum_{j=1}^{N+1}{\left[ \int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x} \right]u_j}=\int_a^b{f\phi _{i}\mathrm{d}x}
,\quad i=1,\cdots,N+1
\end{aligned}
$$
所以
$$
A_{ij}u_j=b_i
,\quad i,j=1,\cdots,N+1
$$
其中
$$
\begin{gather}
A_{ij}=\int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x}
\nonumber
\\
b_i=\int_a^b{f\phi _{i}\mathrm{d}x}
\nonumber
\end{gather}
$$
### Assembly
$$
\begin{gather}
A_{ij}=\int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x}=\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{c\psi _{nj}^{\prime}\psi _{ni}^{\prime}\mathrm{d}x}},i,j=1,\cdots ,N+1
\nonumber
\\
b_i=\int_a^b{f\phi _i\mathrm{d}x}=\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{f\psi _{ni}\mathrm{d}x}},i,j=1,\cdots ,N+1
\nonumber
\end{gather}
$$
## 1.3 Boundary treatment

对于 Dirichlet 边界条件，相当于我们已经知道节点上的精确值，所以那些测试什么的都不要去做了，只需要将组装矩阵的相应部分改成1，其余置零即可，也就是形成一个这样的方程
$$
u_j=g_b
$$
## 1.4 FE Method

1. 生成信息矩阵（顶点信息和单元信息）: $P, T$
2. 组装矩阵和向量: 只基于信息矩阵处理
3. 处理边界条件
4. 求解线性方程组

指导性编程主要在于组装矩阵的子程序编写，也就是这个公式
$$
A_{ij}=\int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x}=\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{c\psi _{nj}^{\prime}\psi _{ni}^{\prime}\mathrm{d}x}},i,j=1,\cdots ,N+1
$$
```matlab
function A = assemble_matrix_1D(?)

A = sparse(matrix_size(1),matrix_size(2));

for n = 1:number_of_elements
    
    vertices = P(:,T(:,n));
    
    [Gauss_weights,Gauss_nodes] = generate_Gauss_local(vertices,Gauss_type);
    
    for alpha = 1:number_of_local_basis_fun_trial
        
        for beta = 1:number_of_local_basis_fun_test
            
            int_value = Gauss_quad_1D_trial_test(?);
            
            A(Tb_test(beta,n),Tb_trial(alpha,n)) = A(Tb_test(beta,n),Tb_trial(alpha,n))+int_value;
            
        end
        
    end
    
end

end
```
打?表示省略参数输入，因为如前面所说，为了便于教学理解，变量命名都写得比较繁琐，同时为了有限元程序的通用性结构，留下的参数接口非常多，比较“又臭又长”。还有这个MATLAB三重大循环，也完全可以向量化处理。

这里主要关注 `vertices = P(:,T(:,n));` 单元顶点信息的数据结构和操作方法，将网格信息和矩阵组装分离开来。

以及积分的子程序 `Gauss_quad_1D_trial_test()` 为了程序的通用性，这里的计算公式是
$$
r=\int_{x_n}^{x_{n+1}}{c\varphi _{n\alpha}^{\left( r \right)}\psi _{n\beta}^{\left( s \right)}}
$$
```matlab
function int_value = Gauss_quad_1D_trial_test(?)

Gauss_point_number = size(Gauss_weights,2);

int_value = 0;

for k = 1:Gauss_point_number
    
    int_value = int_value + Gauss_weights(k)*feval(coe_fun,Gauss_nodes(k))*FE_basis_local_fun_1D(?)*FE_basis_local_fun_1D(?);
    
end

end
```
也就是说编写好一个函数库，通过输入trail function和test function的类型、导数阶等参数，调用函数库来计算，尽管程序会运行很慢，但是保证了通用性。

## 1.5 General extensions

### "Reference -> Local -> Global" Framework

$$
\begin{gather}
\psi \left( x \right) =\hat{\psi}\left( \hat{x} \right) 
\nonumber
\\
\frac{\mathrm{d}\psi _{nj}\left( x \right)}{\mathrm{d}x}=\frac{\mathrm{d}\hat{\psi}_j\left( \hat{x} \right)}{\mathrm{d}\hat{x}}\frac{\mathrm{d}\hat{x}}{\mathrm{d}x}
\nonumber
\end{gather}
$$

### Neumann/Robin boundary conditions

#### Neumann
$$
\begin{gather}
-\frac{\mathrm{d}}{\mathrm{d}x}\left( c\left( x \right) \frac{\mathrm{d}u\left( x \right)}{\mathrm{d}x} \right) =f\left( x \right) , a<x<b
\nonumber
\\
u^{\prime}\left( a \right) =r_a, u\left( b \right) =g_b
\nonumber
\end{gather}
$$
Hence
$$
\int_a^b{cu^{\prime}v^{\prime}\mathrm{d}x=\int_a^b{fv\mathrm{d}x}}-r_ac\left( a \right) v\left( a \right) 
$$
只需在向量阵的合适位置补上 $-r_ac(a)$ 即可，这里是第一项，因为 $v(a)$ 只在 $a$ 处有值，且值为1

#### Robin
$$
\begin{gather}
-\frac{\mathrm{d}}{\mathrm{d}x}\left( c\left( x \right) \frac{\mathrm{d}u\left( x \right)}{\mathrm{d}x} \right) =f\left( x \right) ,a<x<b
\nonumber
\\
u\left( a \right) =g_a,u^{\prime}\left( b \right) +q_bu\left( b \right) =p_b
\nonumber
\end{gather}
$$
Hence
$$
q_bc\left( b \right) u\left( b \right) v\left( b \right) +\int_a^b{cu^{\prime}v^{\prime}\mathrm{d}x}=\int_a^b{fv\mathrm{d}x}+p_bc\left( b \right) v\left( b \right) 
$$
只需在组装矩阵的最后一行最后一列补上 $q_bc(b)$ ，因为 $v(b)$ 表示应在最后一行，$u(b)$ 表示应在最后一列。同时，在向量阵的最后一项补上 $p_bc(b)$ 即可。

### More measurements for errors

- $L^{\infty}$ norm error: $\left\| u-u_h \right\| _{\infty}=\underset{x\in I}{\mathrm{sup}}\left| u\left( x \right)-u_h\left( x \right) \right|\,$
- $L^2$ norm error: $\left\| u-u_h \right\| _0=\sqrt{\int_I{\left( u-u_h \right) ^2\mathrm{d}x}}$
- $H^1$ semi-norm error: $\left| u-u_h \right|_1=\sqrt{\int_I{\left( u^{\prime}-u_{h}^{\prime} \right) ^2\mathrm{d}x}}$

By using $u_h=\sum_{j=1}^{N_b}{u_j\phi _j}$ 
$$
\begin{aligned}
\left\| u-u_h \right\| _{\infty}&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -u_h\left( x \right) \right|
\\
&=\,\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -\sum_{j=1}^{N_b}{u_j\phi _j\left( x \right)} \right|
\\
&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}\left( x \right)} \right|
\end{aligned}
$$
$$
\begin{aligned}
\left\| u-u_h \right\| _0&=\sqrt{\int_I{\left( u-u_h \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-u_h \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{j=1}^{N_b}{u_j\phi _j} \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}} \right) ^2\mathrm{d}x}}}
\end{aligned}
$$
$$
\begin{aligned}
\left| u-u_h \right|_1&=\sqrt{\int_I{\left( u^{\prime}-u_{h}^{\prime} \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u^{\prime}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}^{\prime}} \right) ^2\mathrm{d}x}}}
\end{aligned}
$$
可以发现，误差计算的关键在于计算下列公式的子程序
$$
w_{n,s}=\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}^{\left( s \right)}}
$$
由于基函数库已经写好了，就只需要理解如何从已经计算好的有限元解 $u_h$ 中取出单元上的有限元局部解向量 $u_{T_b\left( k,n \right)}$

## 作业思考

### HW2
改二阶单元主要修改了以下部分：
1. 在库函数中添加了二阶函数，共3个基函数以及各阶导数
2. 修改了Pb，Tb，也就是有限元顶点信息和有限元单元信息，一定要注意区分有限元单元和网格单元这两个不同的概念，尽管一阶单元它们是一样的
3. 最后更改了边界条件处理中犯的错误。A矩阵的行列是与全局有限元节点有关，而不是全局网格节点，HW1写一阶单元没有注意区分

### HW3
写HW3的时候突然想到一个问题，Tb矩阵的排列必须要和基函数的顺序是一致的，也就是说以二阶单元为例，Tb(3,k)是第k个单元的中点，那么基函数index为3的就是中点对应的节点基函数。

同时，在一维单元中没有太多在意边界条件的流程化处理，之后二维单元中再来处理。

# Chapter 2: 2D/3D Finite Element Spaces

## 2.1 2D uniform Mesh

求解区域
$$
\Omega =\left[ left,right \right] \times \left[ bottom,top \right] 
$$
划分网格后的步长
$$
h=\left[ h_1,h_2 \right] =\left[ \frac{right-left}{N_1},\frac{top-bottom}{N_2} \right] 
$$
因此共有 $N=2N_1N_2$ 个矩形网格以及 $N_m=(N_1+1)(N_2+1)$ 个节点

## 2.2 Triangular elements

首先定义仿射变换
$$
J=\left( \begin{matrix}
	x_2-x_1&		x_3-x_1\\
	y_2-y_1&		y_3-y_1\\
\end{matrix} \right) 
$$
$$
\begin{aligned}
\hat{x}&=\frac{\left( y_3-y_1 \right) \left( x-x_1 \right) -\left( x_3-x_1 \right) \left( y-y_1 \right)}{\left| J \right|}
\\
\hat{y}&=\frac{-\left( y_2-y_1 \right) \left( x-x_1 \right) +\left( x_2-x_1 \right) \left( y-y_1 \right)}{\left| J \right|}
\end{aligned}
$$
这样可以使用参考基函数
线性基函数
$$
\begin{aligned}
\hat{\psi}_1\left( \hat{x},\hat{y} \right) &=-\hat{x}-\hat{y}+1
\\
\hat{\psi}_2\left( \hat{x},\hat{y} \right) &=\hat{x}
\\
\hat{\psi}_3\left( \hat{x},\hat{y} \right) &=\hat{y}
\end{aligned}
$$
二阶基函数
$$
\begin{aligned}
\hat{\psi}_1\left( \hat{x},\hat{y} \right) &=2\hat{x}^2+2\hat{y}^2+4\hat{x}\hat{y}-3\hat{x}-3\hat{y}+1
\\
\hat{\psi}_2\left( \hat{x},\hat{y} \right) &=2\hat{x}^2-\hat{x}
\\
\hat{\psi}_3\left( \hat{x},\hat{y} \right) &=2\hat{y}^2-\hat{y}
\\
\hat{\psi}_4\left( \hat{x},\hat{y} \right) &=-4\hat{x}^2-4\hat{x}\hat{y}+4\hat{x}
\\
\hat{\psi}_5\left( \hat{x},\hat{y} \right) &=4\hat{x}\hat{y}
\\
\hat{\psi}_6\left( \hat{x},\hat{y} \right) &=-4\hat{y}^2-4\hat{x}\hat{y}+4\hat{y}
\end{aligned}
$$
采用链式求导法则
$$
\begin{aligned}
\frac{\partial \psi _i}{\partial x}=&\frac{\partial \hat{\psi}_i}{\partial \hat{x}}\frac{y_3-y_1}{\left| J \right|}+\frac{\partial \hat{\psi}_i}{\partial \hat{y}}\frac{y_1-y_2}{\left| J \right|}
\\
\frac{\partial \psi _i}{\partial y}=&\frac{\partial \hat{\psi}_i}{\partial \hat{x}}\frac{x_1-x_3}{\left| J \right|}+\frac{\partial \hat{\psi}_i}{\partial \hat{y}}\frac{x_2-x_1}{\left| J \right|}
\\
\frac{\partial ^2\psi _i}{\partial x^2}=&\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}^2}\frac{\left( y_3-y_1 \right) ^2}{\left| J \right|^2}+2\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}\partial \hat{y}}\frac{\left( y_3-y_1 \right) \left( y_1-y_2 \right)}{\left| J \right|^2}+\frac{\partial ^2\hat{\psi}_i}{\partial \hat{y}^2}\frac{\left( y_1-y_2 \right) ^2}{\left| J \right|^2}
\\
\frac{\partial ^2\psi _i}{\partial y^2}=&\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}^2}\frac{\left( x_1-x_3 \right) ^2}{\left| J \right|^2}+2\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}\partial \hat{y}}\frac{\left( x_1-x_3 \right) \left( x_2-x_1 \right)}{\left| J \right|^2}+\frac{\partial ^2\hat{\psi}_i}{\partial \hat{y}^2}\frac{\left( x_2-x_1 \right) ^2}{\left| J \right|^2}
\\
\frac{\partial ^2\psi _i}{\partial x\partial y}=&\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}^2}\frac{\left( x_1-x_3 \right) \left( y_3-y_1 \right)}{\left| J \right|^2}+\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}\partial \hat{y}}\frac{\left( x_1-x_3 \right) \left( y_1-y_2 \right)}{\left| J \right|^2}
\\
&+\frac{\partial ^2\hat{\psi}_i}{\partial \hat{x}\partial \hat{y}}\frac{\left( x_2-x_1 \right) \left( y_3-y_1 \right)}{\left| J \right|^2}+\frac{\partial ^2\hat{\psi}_i}{\partial \hat{y}^2}\frac{\left( x_2-x_1 \right) \left( y_1-y_2 \right)}{\left| J \right|^2}
\end{aligned}
$$

## 作业思考

网格的索引和有限元的索引可以是不重合的，只需要在使用的时候使用不同的标记即可，在通用程序的编写中，`P, T, Pb, Tb` 都生成了
``` matlab
P(:,T(:,n));
Pb(:,Tb(:,n));
```

又看了一遍之前的程序，P_basis在组装矩阵的时候是不需要用到的，只需要T_basis的组装顺序，这是因为在求解基函数的时候已经用了中点的位置信息。后续可能使用边界节点信息的时候会需要用到P_basis

# Chapter 3: Finite Elements for 2D second order elliptic equation

## 3.1 Weak/Galerkin formulation

考虑二阶椭圆问题
$$
\begin{aligned}
-\nabla \cdot \left( c\nabla u \right) &=f, \ \mathrm{in} \ \Omega 
\\
u&=g,\ \mathrm{on} \  \partial \Omega 
\end{aligned}
$$
其中梯度和散度为
$$
\begin{aligned}
\nabla u &= \left( u_x,u_y \right) 
\\
\nabla \cdot \vec{v} &=\frac{\partial v_1}{\partial x}+\frac{\partial v_2}{\partial y}
\end{aligned}
$$
两边同时乘以测试函数
$$
-\int_{\Omega}{\nabla \cdot \left( c\nabla u \right) v \ \mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv \ \mathrm{d}x\mathrm{d}y}
$$
使用格林公式
$$
\int_{\Omega}{\nabla \cdot \left( c\nabla u \right) v \ \mathrm{d}x\mathrm{d}y}=\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v \ \mathrm{d}s}-\int_{\Omega}{c\nabla u\cdot \nabla v \ \mathrm{d}x\mathrm{d}y}
$$
可以得到
$$
\int_{\Omega}{c\nabla u\cdot \nabla v \ \mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v \ \mathrm{d}s}=\int_{\Omega}{fv \ \mathrm{d}x\mathrm{d}y}
$$
由于Dirichlet条件给定了边界条件，所以取测试函数在边界处为0，因此有
$$
\int_{\Omega}{c\nabla u\cdot \nabla v \ \mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv \ \mathrm{d}x\mathrm{d}y}
$$
和第一章类似，使用一个有限维空间近似无穷维空间
$$
\int_{\Omega}{c\nabla u_h \cdot \nabla v_h \ \mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv_h \ \mathrm{d}x\mathrm{d}y}
$$
## 3.2 FE discretization

取
$$
u_h=\sum_{j=1}^{N_b}{u_j\phi _j}
$$
代入上式
$$
\int_{\Omega}{c\nabla \left( \sum_{j=1}^{N_b}{u_j\phi _j} \right) \cdot \nabla \phi _i\ \mathrm{d}x\mathrm{d}y}=\int_{\Omega}{f\phi _i\ \mathrm{d}x\mathrm{d}y}
$$
整理可得
$$
\sum_{j=1}^{N_b}{u_j}\left[ \int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\ \mathrm{d}x\mathrm{d}y} \right] =\int_{\Omega}{f\phi _i\ \mathrm{d}x\mathrm{d}y}
$$
和第一章类似有
$$
\begin{aligned}
A_{ij}&=\int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\mathrm{d}x\mathrm{d}y}
\\
b_i&=\int_{\Omega}{f\phi _i\mathrm{d}x\mathrm{d}y}
\end{aligned}
$$
最后求解线性代数方程组
$$
A_{ij}u_j=b_i
$$
单元刚度阵组装总刚度阵
$$
A_{ij}=\sum_{n=1}^N{\int_{En}{c\nabla \psi _{n\alpha}\cdot \nabla \psi _{n\beta}\ \mathrm{d}x\mathrm{d}y}}\,\,\left( \alpha ,\beta =1,\cdots ,N_{lb} \right) 
$$
注意到
$$
\int_{En}{c\nabla \psi _{n\alpha}\cdot \nabla \psi _{n\beta}\mathrm{d}x\mathrm{d}y}=\int_{En}{c\frac{\partial \psi _{n\alpha}}{\partial x}\frac{\partial \psi _{n\beta}}{\partial x}\mathrm{d}x\mathrm{d}y}+\int_{En}{c\frac{\partial \psi _{n\alpha}}{\partial y}\frac{\partial \psi _{n\beta}}{\partial y}\mathrm{d}x\mathrm{d}y}
$$
所以更一般的积分形式为
$$
\int_{En}{c\frac{\partial ^{r+s}\psi _{n\alpha}}{\partial x^r\partial y^s}\frac{\partial ^{p+q}\psi _{n\beta}}{\partial x^p\partial y^q}\mathrm{d}x\mathrm{d}y}
$$
## 3.3 Dirichlet boundary condition

和第一章类似

## 3.4 Measurements for errors

- $L^{\infty}$ norm error: $\left\| u-u_h \right\| _{\infty}=\underset{\left( x,y \right) \in {\Omega}}{\mathrm{sup}}\left| u\left( x,y \right)-u_h\left( x,y \right) \right|\,$
- $L^2$ norm error: $\left\| u-u_h \right\| _0=\sqrt{\int_{\Omega}{\left( u-u_h \right) ^2\mathrm{d}x}}$
- $H^1$ semi-norm error: $\left| u-u_h \right|_1=\sqrt{\int_{\Omega}{\left( \frac{\partial \left( u-u_h \right)}{\partial x} \right) ^2\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\left( \frac{\partial \left( u-u_h \right)}{\partial y} \right) ^2\mathrm{d}x\mathrm{d}y}}$

具体计算公式如下
$$
\begin{aligned}
	\left\| u-u_h \right\| _{\infty}&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x,y \right) -u_h\left( x,y \right) \right|\\
	&=\,\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x,y \right) -\sum_{j=1}^{N_b}{u_j\phi _j\left( x,y \right)} \right|\\
	&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x,y \right) -\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}\left( x,y \right)} \right|\\
\end{aligned}
$$
$$
\begin{aligned}
\left\| u-u_h \right\| _0&=\sqrt{\int_I{\left( u-u_h \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-u_h \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{j=1}^{N_b}{u_j\phi _j} \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}} \right) ^2\mathrm{d}x}}}
\end{aligned}
$$
$$
\begin{aligned}
\left| u-u_h \right|_{1}^{2}=&\sum_{n=1}^N{\int_{E_n}{\left( \frac{\partial u}{\partial x}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\frac{\partial \psi _{nk}}{\partial x}} \right) ^2\mathrm{d}x\mathrm{d}y}}
\\
&+\sum_{n=1}^N{\int_{E_n}{\left( \frac{\partial u}{\partial y}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\frac{\partial \psi _{nk}}{\partial y}} \right) ^2\mathrm{d}x\mathrm{d}y}}
\end{aligned}
$$

## 3.5 More Discussion

### 3.5.1 Neumann boundary conditions

$$
\begin{aligned}
\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}&=\int_{\Gamma _N}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}+\int_{\partial \Omega /\Gamma _N}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}
\\
&=\int_{\Gamma _N}{cpv\mathrm{d}s}
\end{aligned}
$$
离散后的有限元格式为
$$
\int_{\Omega}{c\nabla u_h\cdot \nabla v_h\,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv_h\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _N}{cpv_h\mathrm{d}s}
$$
整理后有
$$
\sum_{j=1}^{N_b}{u_j}\left[ \int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] =\int_{\Omega}{f\phi _i\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _N}{cp\phi _i\mathrm{d}s}
$$
所以更新右端载荷矩阵即可
$$
v_i=\int_{\Gamma _N}{cp\phi _i\mathrm{d}s}
\\
=\sum_{1\le k\le nbe}{\int_{e_k}{cp\psi _{n_k\beta}\mathrm{d}s}}\,\,\left( \beta =1,\cdots ,N_{lb} \right) 
$$
### 3.5.2 Robin boundary conditions
边界条件为
$$
\nabla u\cdot \vec{n}+ru=q
$$
带入原式有
$$
\begin{aligned}
	\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}&=\int_{\Gamma _R}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}+\int_{\partial \Omega /\Gamma _R}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}\\
	&=\int_{\Gamma _R}{c\left( q-ru \right) v\mathrm{d}s}\\
	&=\int_{\Gamma _R}{cqv\mathrm{d}s}-\int_{\Gamma _R}{cruv\mathrm{d}s}\\
\end{aligned}
$$
离散后的有限元格式为
$$
\int_{\Omega}{c\nabla u_h\cdot \nabla v_h\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _R}{cru_hv_h\mathrm{d}s}=\int_{\Omega}{fv_h\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _R}{cpv_h\mathrm{d}s}
$$
整理后有
$$
\sum_{j=1}^{N_b}{u_j}\left[ \int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] +\sum_{j=1}^{N_b}{u_j}\left[ \int_{\Gamma _R}{cr\phi _j\phi _i\,\,\mathrm{d}s} \right] =\int_{\Omega}{f\phi _i\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _N}{cp\phi _i\mathrm{d}s}
$$
刚度阵和载荷阵都需要更新！

## 作业思考

二维问题修改的最大工作量在上一章的作业，也就是网格信息、有限元单元信息、边界边信息和边界节点信息的生成。其他差不多只需要将有x的地方改成x和y就行了，也就是这个公式的计算子程序编写
$$
\int_{En}{c\frac{\partial ^{r+s}\psi _{n\alpha}}{\partial x^r\partial y^s}\frac{\partial ^{p+q}\psi _{n\beta}}{\partial x^p\partial y^q}\mathrm{d}x\mathrm{d}y}
$$

这里需要强调一点，边界边信息和边界节点信息分属于不同类别，边界边为网格概念，边界节点为有限元概念。我在上一章的作业中混淆了这两个概念，刚好后面何老师在答疑课中也讲了这个问题，具体链接如下，时间大约在36:39[有限元基础编程课程答疑-何晓明-2021-1-07](https://www.bilibili.com/video/BV1ET4y1N77h/?share_source=copy_web&vd_source=4846b273dd993443973d15943a9b6546)

### `boundaryNodes` and `boundaryEdges`

- `boundaryNodes` 存边界信息，所有边界上的有限元节点
- `boundaryEdges` 存边界信息，这条边所在的有限元单元索引，以及边上的两个顶点（网格点）用来确定边的位置
> 感觉这里可以设计更好的数据结构，在复习重构的时候考虑一下这个问题

# Chapter 4: Finite Elements for 2D second order parabolic and hyperbolic equation

## 4.1 Weak formulation /4.2 Semi-discretization

Cosider the 2D second order parabolic equation
$$
\begin{aligned}
&u_t-\nabla \cdot \left( c\nabla u \right) =f, \  \mathrm{in} \ \Omega \times \left[ 0,T \right] 
\\
&u=g,\ \mathrm{on}\ \partial \Omega \times \left[ 0,T \right] 
\\
&u=u_0,\ \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega 
\end{aligned}
$$
首先右乘测试函数并积分
$$
\int_{\Omega}{u_tv\,\,\mathrm{d}x\mathrm{d}y}-\int_{\Omega}{\nabla \cdot \left( c\nabla u \right) v\,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv\,\,\mathrm{d}x\mathrm{d}y}
$$
使用格林公式整理可得
$$
\int_{\Omega}{u_tv\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{c\nabla u\cdot \nabla v\,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v\,\,\mathrm{d}s}=\int_{\Omega}{fv\,\,\mathrm{d}x\mathrm{d}y}
$$
考虑Dirichlet边界条件，有限元离散
$$
\int_{\Omega}{u_{h_t}v_h\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{c\nabla u_h\cdot \nabla v_h\,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{fv_h\,\,\mathrm{d}x\mathrm{d}y}
$$
离散格式为
$$
u_h\left( x,y,t \right) =\sum_{j=1}^{N_b}{u_j\left( t \right) \phi _j\left( x,y \right)}
$$
代入整理
$$
\sum_{j=1}^{N_b}{\dot{u}_j\left( t \right)}\left[ \int_{\Omega}{\phi _j\phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] +\sum_{j=1}^{N_b}{u_j\left( t \right)}\left[ \int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] =\int_{\Omega}{f\phi _i\,\,\mathrm{d}x\mathrm{d}y}
$$
最后简化为ODE
$$
M_{ij}\dot{u}_j+A_{ij}u_j=b_i
$$
## 4.3 Full discretization

这里只抄写最后的差分格式 $\theta$-scheme
$$
\frac{y_{j+1}-y_j}{h}=\theta f\left( t_{j+1},y_{j+1} \right) +\left( 1-\theta \right) f\left( t_j,y_j \right) 
$$
将之前得到的ODE写成差分格式
$$
M\frac{u_{m+1}-u_m}{\Delta t}+\theta A\left( t_{m+1} \right) u_{m+1}+\left( 1-\theta \right) A\left( t_m \right) u_m=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) 
$$
最后写成线性方程组，每一个时间步求解一次线性方程组就可以了
$$
\tilde{A}_{m+1}u_{m+1}=\tilde{b}_{m+1}
$$
其中
$$
\begin{aligned}
\tilde{A}_{m+1}&=\frac{M}{\Delta t}+\theta A\left( t_{m+1} \right) 
\\
\tilde{b}_{m+1}&=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) +\left( \frac{M}{\Delta t}-\left( 1-\theta \right) A\left( t_m \right) \right) u_m
\end{aligned}
$$
简单起见，考虑LTI系统
$$
\begin{aligned}
\tilde{A}&=\frac{M}{\Delta t}+\theta A
\\
\tilde{b}_{m+1}&=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) +\left[ \frac{M}{\Delta t}-\left( 1-\theta \right) A \right] u_m
\end{aligned}
$$

# Chapter 5: Finite elements for 2D steady linear elasticity equation

## 5.1 Weak/Galerkin formulation

Consider the 2D linear elasticity equation:
$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,   \mathrm{in}\ \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,   \mathrm{on}\ \partial \Omega\\
\end{cases}
$$
应力张量定义为
$$
\sigma \left( \boldsymbol{u} \right) =\left( \begin{matrix}
	\sigma _{11}\left( \boldsymbol{u} \right)&		\sigma _{12}\left( \boldsymbol{u} \right)\\
	\sigma _{21}\left( \boldsymbol{u} \right)&		\sigma _{22}\left( \boldsymbol{u} \right)\\
\end{matrix} \right) ,\quad \sigma _{ij}\left( \boldsymbol{u} \right) =\lambda \left( \nabla \cdot \boldsymbol{u} \right) \delta _{ij}+2\mu \epsilon _{ij}\left( \boldsymbol{u} \right) 
$$
应变张量定义为
$$
\epsilon =\left( \begin{matrix}
	\epsilon _{11}&		\epsilon _{12}\\
	\epsilon _{21}&		\epsilon _{22}\\
\end{matrix} \right) ,\quad \epsilon _{ij}\left( \boldsymbol{u} \right) =\frac{1}{2}\left( \frac{\partial u_i}{\partial x_j}+\frac{\partial u_j}{\partial x_i} \right) 
$$
整理可得
$$
\sigma \left( \boldsymbol{u} \right) =\left( \begin{matrix}
	\lambda \frac{\partial u_1}{\partial x_1}+2\mu \frac{\partial u_1}{\partial x_1}+\lambda \frac{\partial u_2}{\partial x_2}&		\mu \frac{\partial u_1}{\partial x_2}+\mu \frac{\partial u_2}{\partial x_1}\\
	\mu \frac{\partial u_1}{\partial x_2}+\mu \frac{\partial u_2}{\partial x_1}&		\lambda \frac{\partial u_1}{\partial x_1}+2\mu \frac{\partial u_2}{\partial x_2}+\lambda \frac{\partial u_2}{\partial x_2}\\
\end{matrix} \right) 
$$
等式两边同时点乘测试函数并积分
$$
-\int_{\Omega}{\left( \nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) \right) \cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
又因为
$$
\int_{\Omega}{\left( \nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) \right) \cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}=\int_{\partial \Omega}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}-\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right): \nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
整理可得
$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right): \nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}-\int_{\partial \Omega}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
考虑Dirichlet边界条件，写成标量形式
$$
\begin{aligned}
&\int_{\Omega}{ ( \lambda \frac{\partial u_1}{\partial x_1}\frac{\partial v_1}{\partial x_1}+2\mu \frac{\partial u_1}{\partial x_1}\frac{\partial v_1}{\partial x_1}+\lambda \frac{\partial u_2}{\partial x_2}\frac{\partial v_1}{\partial x_1}\,\,}
\\
&+\mu \frac{\partial u_1}{\partial x_2}\frac{\partial v_1}{\partial x_2}+\mu \frac{\partial u_2}{\partial x_1}\frac{\partial v_1}{\partial x_2}+\mu \frac{\partial u_1}{\partial x_2}\frac{\partial v_2}{\partial x_1}+\mu \frac{\partial u_2}{\partial x_1}\frac{\partial v_2}{\partial x_1}
\\
&+\lambda \frac{\partial u_1}{\partial x_1}\frac{\partial v_2}{\partial x_2}+\lambda \frac{\partial u_2}{\partial x_2}\frac{\partial v_2}{\partial x_2}+2\mu \frac{\partial u_2}{\partial x_2}\frac{\partial v_2}{\partial x_2}  )\,\,\mathrm{d}x_1\mathrm{d}x_2
\\
=&\int_{\Omega}{\left( f_1v_1+f_2v_2 \right) \,\, \mathrm{d}x_1\mathrm{d}x_2}
\end{aligned}
$$
有限元空间逼近
$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u}_h \right) :\nabla \boldsymbol{v}_h\,\,\mathrm{d}x_1\mathrm{d}x_2}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}_h\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
对于解有
$$
u_{1h}=\sum_{j=1}^{N_b}{u_{1j}\phi _j},   u_{2h}=\sum_{j=1}^{N_b}{u_{2j}\phi _j}
$$
可以看到，未知数为原来的两倍，因此需要更多的测试函数。首先取 $\boldsymbol{v}_h=\left( \phi _i,0 \right) ^T$ 
$$
\begin{aligned}
&\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}+2\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&=\int_{\Omega}{f_1\phi _i\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{aligned}
$$
再取 $\boldsymbol{v}_h=\left( 0,\phi _i \right) ^T$ 
$$
\begin{aligned}
&\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}+2\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&=\int_{\Omega}{f_2\phi _i\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{aligned}
$$
排列成矩阵形式
$$
\begin{aligned}
A_1=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2},\quad    A_2=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_3=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2},  \quad   A_4=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_5=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2},    \quad A_6=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_7=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2},    \quad A_8=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{aligned}
$$
写为
$$
\left( \begin{matrix}
	A_1+2A_2+A_3&		A_4+A_5\\
	A_6+A_7&		A_8+2A_3+A_2\\
\end{matrix} \right) \left( \begin{array}{c}
	u_1\\
	u_2\\
\end{array} \right) =\left( \begin{array}{c}
	b_1\\
	b_2\\
\end{array} \right) 
$$
也就是说，只需要多次调用矩阵组装的子程序，最后排成大矩阵求解即可

## 5.2 Measurements for errors

$L^{\infty}$ norm error:
$$
\begin{aligned}
\left\| \boldsymbol{u}-\boldsymbol{u}_h \right\| _{\infty}
&=\max \left( \left\| u_1-u_{1h} \right\| _{\infty},\left\| u_2-u_{2h} \right\| _{\infty} \right) 
\\
\left\| u_1-u_{1h} \right\| _{\infty}
&=\underset{\Omega}{\mathrm{sup}}\left| u_1-u_{1h} \right|
\\
\left\| u_2-u_{2h} \right\| _{\infty}
&=\underset{\Omega}{\mathrm{sup}}\left| u_2-u_{2h} \right|
\end{aligned}
$$


$L^2$ norm error:
$$
\begin{aligned}
\left\| \boldsymbol{u}-\boldsymbol{u}_h \right\| _0&=\sqrt{\left\| u_1-u_{1h} \right\| _{0}^{2}+\left\| u_2-u_{2h} \right\| _{0}^{2}}
\\
\left\| u_1-u_{1h} \right\| _0&=\sqrt{\int_{\Omega}{\left( u_1-u_{1h} \right) ^2\mathrm{d}x_1\mathrm{d}x_2}}
\\
\left\| u_2-u_{2h} \right\| _0&=\sqrt{\int_{\Omega}{\left( u_2-u_{2h} \right) ^2\mathrm{d}x_1\mathrm{d}x_2}}
\end{aligned}
$$

$H^1$ semi-norm error:
$$
\begin{aligned}
\left| \boldsymbol{u}-\boldsymbol{u}_h \right|_1&=\sqrt{\left| u_1-u_{1h} \right|_{1}^{2}+\left| u_2-u_{2h} \right|_{1}^{2}}
\\
\left| u_1-u_{1h} \right|_1&=\sqrt{\int_{\Omega}{\left( \frac{\partial \left( u_1-u_{1h} \right)}{\partial x_1} \right) ^2+\left( \frac{\partial \left( u_1-u_{1h} \right)}{\partial x_2} \right) ^2\mathrm{d}x_1\mathrm{d}x_2}}
\\
\left| u_2-u_{2h} \right|_1&=\sqrt{\int_{\Omega}{\left( \frac{\partial \left( u_2-u_{2h} \right)}{\partial x_1} \right) ^2+\left( \frac{\partial \left( u_2-u_{2h} \right)}{\partial x_2} \right) ^2\mathrm{d}x_1\mathrm{d}x_2}}
\end{aligned}
$$

## 5.3 Stress boundary condition

$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,   \mathrm{in}\ \Omega\\
	\sigma (\boldsymbol{u}) \boldsymbol{n}=\boldsymbol{p}\,\, \mathrm{on}\ \partial\Omega\\
\end{cases}
$$
根据之前得到的公式
$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right): \nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}-\int_{\partial \Omega}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
可以直接将边界条件代入
$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right): \nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}-\int_{\partial \Omega}{\boldsymbol{p} \cdot \boldsymbol{v}\,\,\mathrm{d}s}
=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}
$$
注意在这种边界条件下解不唯一，需要加入 Dirichlet 边界条件，因此考虑下面的混合边界条件
$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\sigma (\boldsymbol{u})\boldsymbol{n}=\boldsymbol{p}\, \, \mathrm{on}\ \Gamma _S\subset \partial \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \Gamma _D=\partial \Omega /\Gamma _S\\
\end{cases}
$$
其中 $\boldsymbol{n}=(n_1,n_2)^T$ 是 $\Gamma_S$ 的单位外法线向量

因此边界积分项可以写为
$$
\begin{aligned}
\int_{\partial \Omega}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}&=\int_{\Gamma _S}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}+\int_{\partial \Omega /\Gamma _S}{\left( \mathbf{\sigma }\left( \boldsymbol{u} \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}
\\
&=\int_{\Gamma _S}{\boldsymbol{p}\cdot \boldsymbol{v}\,\,\mathrm{d}s}
\end{aligned}
$$
弱形式为
$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right) :\nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _S}{\boldsymbol{p}\cdot \boldsymbol{v}\,\,\mathrm{d}s}
$$

## 5.3 Robin boundary conditions

$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\sigma (\boldsymbol{u})\boldsymbol{n}+r\boldsymbol{u}=\boldsymbol{q}\,\,\mathrm{on}\ \Gamma _R\subseteq \partial \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \Gamma _D=\partial \Omega /\Gamma _S\\
\end{cases}
$$

$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right) :\nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _R}{r\boldsymbol{u}\cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _R}{\boldsymbol{q}\cdot \boldsymbol{v}\,\,\mathrm{d}s}
$$

## 5.4 Stress boundary condition in normal/tangential directions

$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\boldsymbol{n}^T\sigma (\boldsymbol{u})\boldsymbol{n}=p_n,\  \boldsymbol{\tau }^T\sigma (\boldsymbol{u})\boldsymbol{n}=p_{\tau}\,\,\mathrm{on}\ \Gamma _S\subset \partial \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \Gamma _D=\partial \Omega /\Gamma _S\\
\end{cases}
$$

$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right) :\nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _S}{p_n\left( \boldsymbol{n}^T\boldsymbol{v} \right) \,\,\mathrm{d}s}+\int_{\Gamma _S}{p_{\tau}\left( \boldsymbol{\tau }^T\boldsymbol{v} \right) \,\,\mathrm{d}s}
$$

## 5.5 Robin boundary condition in normal/tangential directions
$$
\begin{cases}
	-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\boldsymbol{n}^T\sigma (\boldsymbol{u})\boldsymbol{n}+r\boldsymbol{n}^T\boldsymbol{u}=q_n,\ \boldsymbol{\tau }^T\sigma (\boldsymbol{u})\boldsymbol{n}+r\boldsymbol{\tau }^T\boldsymbol{u}=q_{\tau}\,\,\mathrm{on}\ \Gamma _R\subseteq \partial \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \Gamma _D=\partial \Omega /\Gamma _R\\
\end{cases}
$$

$$
\int_{\Omega}{\mathbf{\sigma }\left( \boldsymbol{u} \right) :\nabla \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _R}{\left( r\boldsymbol{n}^T\boldsymbol{u} \right) \left( \boldsymbol{n}^T\boldsymbol{v} \right) \,\,\mathrm{d}s}+\int_{\Gamma _R}{\left( r\boldsymbol{\tau }^T\boldsymbol{u} \right) \left( \boldsymbol{\tau }^T\boldsymbol{v} \right) \,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Gamma _R}{q_n\left( \boldsymbol{n}^T\boldsymbol{v} \right) \,\,\mathrm{d}s}+\int_{\Gamma _R}{q_{\tau}\left( \boldsymbol{\tau }^T\boldsymbol{v} \right) \,\,\mathrm{d}s}
$$

# 学习小结

在最简单的一维问题中打好基础，理解好有限元通用结构后，后续的通用程序大部分部分基本不变，每次只需要做少量的修改。然后逐渐向不同的有限元、不同的边界条件、二维三维、非稳态、非线性情况去拓展自己写出来的程序包。

目前完成的部分算是《有限元基础编程I》，差不多是正常教学一学期的内容，后续的内容《有限元基础编程II》录播中也有，之后打算进一步学习。

同时，前半段的学习算是手把手带上路，包括问题描述、公式推导、伪代码课件中都有详细的讲述，指导性编程也比较详细，之后的学习应尝试自己独立走完整个流程，完成课件中提到的一些小课题等内容，自行设计和完成相应的算例测试。

# Chapter 6: Finite elements for 2D steady Stokes equation

## 6.1 Weak/Galerkin formulation

$$
\begin{cases}
	-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \partial \Omega\\
\end{cases}
$$

where
$$
\boldsymbol{u}\left( x,y \right) =\left( u_1,u_2 \right) ^T,\boldsymbol{g}\left( x,y \right) =\left( g_1,g_2 \right) ^T,\boldsymbol{f}\left( x,y \right) =\left( f_1,f_2 \right) ^T
$$

the stress tensor $\mathbb{T} \left( \boldsymbol{u},p \right)$ is defined as
$$
\mathbb{T} \left( \boldsymbol{u},p \right) =2\nu \mathbb{D} \left( \boldsymbol{u} \right) -p\mathbb{I} 
$$

where $\nu$ is the viscosity and the deformation tensor 
$$
\mathbb{D} \left( \boldsymbol{u} \right) =\frac{1}{2}\left( \nabla \boldsymbol{u}+\left( \nabla \boldsymbol{u} \right) ^T \right) 
$$

in more details, the deformation tensor can be written as
$$
\mathbb{D} \left( \boldsymbol{u} \right) =\left[ \begin{matrix}
	\frac{\partial u_1}{\partial x}&		\frac{1}{2}\left( \frac{\partial u_1}{\partial y}+\frac{\partial u_2}{\partial x} \right)\\
	\frac{1}{2}\left( \frac{\partial u_1}{\partial y}+\frac{\partial u_2}{\partial x} \right)&		\frac{\partial u_2}{\partial y}\\
\end{matrix} \right] 
$$

hence the stress tensor can be written as 
$$
\mathbb{T} \left( \boldsymbol{u},p \right) =\left[ \begin{matrix}
	2\nu \frac{\partial u_1}{\partial x}-p&		\nu \left( \frac{\partial u_1}{\partial y}+\frac{\partial u_2}{\partial x} \right)\\
	\nu \left( \frac{\partial u_1}{\partial y}+\frac{\partial u_2}{\partial x} \right)&		2\nu \frac{\partial u_1}{\partial x}-p\\
\end{matrix} \right] 
$$

First, take the inner product with a vector function $\boldsymbol{v}\left( x,y \right) =\left( v_1,v_2 \right) ^T$ on both sides of the Stokes equation:
$$
-\int_{\Omega}{\left( \nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) \right) \cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
$$

Second, multiply the divergence free equation by a function $q(x,y)$:
$$
\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}
$$

$\boldsymbol{u}(x,y)$ and $p(x,y)$ are called trail functions, $\boldsymbol{v}(x,y)$ and $q(x,y)$ are called test functions.

Using integration by parts in multi-dimension:
$$
\int_{\Omega}{\left( \nabla \cdot \mathbb{T} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}=\int_{\partial \Omega}{\left( \mathbb{T} \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}-\int_{\Omega}{\mathbb{T} :\nabla \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
$$

where $\boldsymbol{n}=(n_1,n_2)^T$ is the unit outer normal vector of $\partial \Omega$, we obtain
$$
\int_{\Omega}{\mathbb{T} \left( \boldsymbol{u},p \right) :\nabla \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( \mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
$$

Here,
$$
\begin{aligned}
	\mathbf{A}:\mathbf{B}&=\left[ \begin{matrix}
	a_{11}&		a_{12}\\
	a_{21}&		a_{22}\\
\end{matrix} \right] :\left[ \begin{matrix}
	b_{11}&		b_{12}\\
	b_{21}&		b_{22}\\
\end{matrix} \right] 
\\
&=a_{11}b_{11}+a_{12}b_{12}+a_{21}b_{21}+a_{22}b_{22}
\end{aligned}
$$

At last, we obtain:
$$
\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\Omega}{p\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( \mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
\\
-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0
$$

Here we multiply the second equation by $-1$ in order to keep the matrix formulation symmetric later.

$$
\begin{aligned}
&\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x \mathrm{d}y}
\\
=&\int_{\Omega}{\nu \left( 2\frac{\partial u_1}{\partial x}\frac{\partial v_1}{\partial x}+2\frac{\partial u_2}{\partial y}\frac{\partial v_2}{\partial y}+\frac{\partial u_1}{\partial y}\frac{\partial v_1}{\partial y}+\frac{\partial u_1}{\partial y}\frac{\partial v_2}{\partial x}+\frac{\partial u_2}{\partial x}\frac{\partial v_1}{\partial y}+\frac{\partial u_2}{\partial x}\frac{\partial v_2}{\partial x} \right) \,\mathrm{d}x\mathrm{d}y}
\end{aligned}
$$

## 6.2 Matrix formulation
Define
$$
\mathbf{A}_1=\left[ \int_{\Omega}{\nu \frac{\partial \phi _j}{\partial x}\frac{\partial \phi _i}{\partial x}\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}, \mathbf{A}_2=\left[ \int_{\Omega}{\nu \frac{\partial \phi _j}{\partial y}\frac{\partial \phi _i}{\partial y}\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$

$$
\mathbf{A}_3=\left[ \int_{\Omega}{\nu \frac{\partial \phi _j}{\partial x}\frac{\partial \phi _i}{\partial y}\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}, \mathbf{A}_4=\left[ \int_{\Omega}{\nu \frac{\partial \phi _j}{\partial y}\frac{\partial \phi _i}{\partial x}\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$

$$
\mathbf{A}_5=\left[ \int_{\Omega}{-\psi _j\frac{\partial \phi _i}{\partial x}\mathrm{d}x\mathrm{d}y} \right] _{i=1,j=1}^{N_b,N_{bp}}, \mathbf{A}_6=\left[ \int_{\Omega}{-\psi _j\frac{\partial \phi _i}{\partial y}\mathrm{d}x\mathrm{d}y} \right] _{i=1,j=1}^{N_b,N_{bp}}
$$

$$
\mathbf{A}_7=\left[ \int_{\Omega}{-\frac{\partial \phi _j}{\partial x}\psi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1,j=1}^{N_{bp},N_b}, \mathbf{A}_8=\left[ \int_{\Omega}{-\frac{\partial \phi _j}{\partial y}\psi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1,j=1}^{N_{bp},N_b}
$$

Define a zero matrix $\mathbb{O} _1=\left[ 0 \right] _{i=1,j=1}^{N_{bp},N_{bp}}$ whose size is $N_{bp}\times N_{bp}$, then
$$
\mathbf{A}=\left[ \begin{matrix}
	2\mathbf{A}_1+\mathbf{A}_2&		\mathbf{A}_3&		\mathbf{A}_5\\
	\mathbf{A}_4&		2\mathbf{A}_2+\mathbf{A}_1&		\mathbf{A}_6\\
	\mathbf{A}_7&		\mathbf{A}_8&		\mathbb{O} _1\\
\end{matrix} \right] 
$$
$$
\mathbf{A}_4=\mathbf{A}_{3}^{T}, \mathbf{A}_7=\mathbf{A}_{5}^{T}, \mathbf{A}_8=\mathbf{A}_{6}^{T}
$$

Hence the matrix $\mathbf{A}$ is actually symmetric
$$
\mathbf{A}=\left[ \begin{matrix}
	2\mathbf{A}_1+\mathbf{A}_2&		\mathbf{A}_3&		\mathbf{A}_5\\
	\mathbf{A}_{3}^{T}&		2\mathbf{A}_2+\mathbf{A}_1&		\mathbf{A}_6\\
	\mathbf{A}_{5}^{T}&		\mathbf{A}_{6}^{T}&		\mathbb{O} _1\\
\end{matrix} \right] 
$$

## 6.3 Stress boundary condition

$$
\begin{cases}
	-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega\\
	\mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n}=\boldsymbol{p}\,\,\mathrm{on}\ \Gamma _S\subset \partial \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \Gamma _D=\partial \Omega /\Gamma _S\\
\end{cases}
$$

where 
$$
\boldsymbol{u}\left( x,y \right) =\left( u_1,u_2 \right) ^T,\boldsymbol{g}\left( x,y \right) =\left( g_1,g_2 \right) ^T,\boldsymbol{p}\left( x,y \right) =\left( p_1,p_2 \right) ^T,\boldsymbol{f}\left( x,y \right) =\left( f_1,f_2 \right) ^T
$$

The weak formulation
$$
\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\Omega}{p\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Gamma _S}{\boldsymbol{p}\cdot \boldsymbol{v}\,\,\mathrm{d}s}
\\
-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0
$$

# Chapter 7: Finite elements for 2D steady Navier-Stokes equation

## 7.1 Weak/Galerkin formulation

Consider the 2D Navier-Stokes equation:
$$
\begin{cases}
	\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega\\
	\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \partial \Omega\\
\end{cases}
$$

where
$$
\boldsymbol{u}\left( x,y \right) =\left( u_1,u_2 \right) ^T,\boldsymbol{g}\left( x,y \right) =\left( g_1,g_2 \right) ^T,\boldsymbol{f}\left( x,y \right) =\left( f_1,f_2 \right) ^T
$$

The nonlinear advection is defined as
$$
\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}=\left( \begin{array}{c}
	u_1\frac{\partial u_1}{\partial x}+u_2\frac{\partial u_1}{\partial y}\\
	u_1\frac{\partial u_2}{\partial x}+u_2\frac{\partial u_2}{\partial y}\\
\end{array} \right) 
$$

Weak formulation in the vector format:
$$
\begin{aligned}
&\int_{\Omega}{\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}
\\
&-\int_{\Omega}{p\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( \mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y},
\\
&-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0
\end{aligned}
$$

nonlinear advection in the scalar format:
$$
\begin{aligned}
&\int_{\Omega}{\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
\\
&=\int_{\Omega}{\left( u_1\frac{\partial u_1}{\partial x}v_1+u_2\frac{\partial u_1}{\partial y}v_1+u_1\frac{\partial u_2}{\partial x}v_2+u_2\frac{\partial u_2}{\partial y}v_2 \right) \,\,\mathrm{d}x\mathrm{d}y}
\end{aligned}
$$

## 7.2 Newton's iteration

- Initial guess: $\boldsymbol{u}^{\left( 0 \right)}$ and $p^{(0)}$
- Newton's iteration for the weak formulation in the vector
format: for $l=1,2,...,L$
$$
\begin{aligned}
	&\int_{\Omega}{\left( \boldsymbol{u}^{\left( l \right)}\cdot \nabla \right) \boldsymbol{u}^{\left( l-1 \right)}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\left( \boldsymbol{u}^{\left( l-1 \right)}\cdot \nabla \right) \boldsymbol{u}^{\left( l \right)}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
	\\
	&+\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u}^{\left( l \right)} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\Omega}{p^{\left( l \right)}\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}
	\\
	&=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\left( \boldsymbol{u}^{\left( l-1 \right)}\cdot \nabla \right) \boldsymbol{u}^{\left( l-1 \right)}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}
	\\
	&-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u}^{\left( l \right)} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0\\
\end{aligned}
$$

## 7.3 Matrix formulation

$$
AN_1=\left[ \int_{\Omega}{\frac{\partial u_{1h}^{\left( l-1 \right)}}{\partial x}\phi _j\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}, AN_2=\left[ \int_{\Omega}{u_{1h}^{\left( l-1 \right)}\frac{\partial \phi _j}{\partial x}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$
$$
AN_3=\left[ \int_{\Omega}{u_{2h}^{\left( l-1 \right)}\frac{\partial \phi _j}{\partial y}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}, AN_4=\left[ \int_{\Omega}{\frac{\partial u_{1h}^{\left( l-1 \right)}}{\partial y}\phi _j\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$
$$
AN_5=\left[ \int_{\Omega}{\frac{\partial u_{2h}^{\left( l-1 \right)}}{\partial x}\phi _j\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}, AN_6=\left[ \int_{\Omega}{\frac{\partial u_{2h}^{\left( l-1 \right)}}{\partial y}\phi _j\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$

Define a zero matrix $\mathbb{O}_1=\left[ 0 \right] _{i=1,j=1}^{N_{bp},N_{bp}}$ whose size is $N_{bp}\times N_{bp}$ , $\mathbb{O}_2=\left[ 0 \right] _{i=1,j=1}^{N_{b},N_{bp}}$ whose size is $N_{b}\times N_{bp}$ , then
$$
AN=\left[ \begin{matrix}
	AN_1+AN_2+AN_3&		AN_4&		\mathbb{O} _2\\
	AN_5&		AN_6+AN_2+AN_3&		\mathbb{O} _2\\
	\mathbb{O} _{2}^{T}&		\mathbb{O} _{2}^{T}&		\mathbb{O} _1\\
\end{matrix} \right]
$$

Define the vector
$$
\overrightarrow{bN}=\left[ \begin{array}{c}
	\overrightarrow{bN_1}+\overrightarrow{bN_2}\\
	\overrightarrow{bN_3}+\overrightarrow{bN_4}\\
	\overrightarrow{0}\\
\end{array} \right] 
$$
where
$$
\overrightarrow{bN_1}=\left[ \int_{\Omega}{u_{1h}^{\left( l-1 \right)}\frac{\partial u_{1h}^{\left( l-1 \right)}}{\partial x}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1}^{N_b},\overrightarrow{bN_2}=\left[ \int_{\Omega}{u_{2h}^{\left( l-1 \right)}\frac{\partial u_{1h}^{\left( l-1 \right)}}{\partial y}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1}^{N_b}
$$
$$
\overrightarrow{bN_3}=\left[ \int_{\Omega}{u_{1h}^{\left( l-1 \right)}\frac{\partial u_{2h}^{\left( l-1 \right)}}{\partial x}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1}^{N_b},\overrightarrow{bN_4}=\left[ \int_{\Omega}{u_{2h}^{\left( l-1 \right)}\frac{\partial u_{2h}^{\left( l-1 \right)}}{\partial y}\phi _i\mathrm{d}x\mathrm{d}y} \right] _{i=1}^{N_b}
$$

Here the size of the zero vector is $N_{bp}\times 1$ . That is, $\overrightarrow{0}=\left[ 0 \right] _{i=1}^{N_{bp}}$

Define the unknown vector
$$
\vec{X}^{\left( l \right)}=\left[ \begin{array}{c}
	\vec{X}_{1}^{\left( l \right)}\\
	\vec{X}_{2}^{\left( l \right)}\\
	\vec{X}_{3}^{\left( l \right)}\\
\end{array} \right] 
$$

where
$$
\vec{X}_{1}^{\left( l \right)}=\left[ u_{1j}^{\left( l \right)} \right] _{j=1}^{N_b},\ \vec{X}_{2}^{\left( l \right)}=\left[ u_{2j}^{\left( l \right)} \right] _{j=1}^{N_b},\ \vec{X}_{3}^{\left( l \right)}=\left[ p_{j}^{\left( l \right)} \right] _{j=1}^{N_{bp}}
$$

Define
$$
A^{\left( l \right)}=A+AN,\ \vec{b}^{\left( l \right)}=\vec{b}+\overrightarrow{bN}
$$

For step $l(l=1,2,...,L)$ of the Newton's iteration, we obtain the linear algebratic system
$$
A^{\left( l \right)}\vec{X}^{\left( l \right)}=\vec{b}^{\left( l \right)}
$$

# Chapter 8: Finite elements for 2D unsteady Stokes and linear elasticity equations

## 8.1 Weak formulation

Consider the 2D unsteady Stokes equation
$$
\begin{cases}
	\boldsymbol{u}_t-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega \times \left[ 0,T \right]\\
	\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \partial \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{u}_0, p=p_0, \mathrm{at}\ t=0 \ \mathrm{and}\ \mathrm{in}\ \Omega\\
\end{cases}
$$

where $\Omega$ is a 2D domain, $[0,T]$ is the time interval, $\boldsymbol{f}(x,y,t)$ is a given function on $\Omega\times[0,T]$ , $\boldsymbol{g}(x,y,t)$ is a given function on $\partial\Omega\times[0,T]$ , $\boldsymbol{u}_0(x,y)$ and $p_0(x,y)$ are given functions on $\Omega$ at $t=0$ , $\boldsymbol{u}(x,y,t)$ and $p(x,y,t)$ are the unknown functions, and
$$
\boldsymbol{u}\left( x,y,t \right) =\left( u_1,u_2 \right) ^T,\boldsymbol{g}\left( x,y,t \right) =\left( g_1,g_2 \right) ^T,\\
\boldsymbol{f}\left( x,y,t \right) =\left( f_1,f_2 \right) ^T,
\boldsymbol{u}_0\left( x,y \right) =\left( u_{10},u_{20} \right) ^T,
$$

Weak formulation
$$
\begin{aligned}
	&\int_{\Omega}{\boldsymbol{u}_t\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}\\
	&-\int_{\Omega}{p\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( \mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y,}\\
	&-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0\\
\end{aligned}
$$

We also have
$$
\int_{\Omega}{\boldsymbol{u}_t\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}=\int_{\Omega}{\frac{\partial u_1}{\partial t}v_1\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\frac{\partial u_2}{\partial t}v_2\,\,\mathrm{d}x\mathrm{d}y}
$$

## 8.2 Matrix formulation

Define the basic mass matrix
$$
M_e=\left[ m_{ij} \right] _{i,j=1}^{N_b}=\left[ \int_{\Omega}{\phi _j\phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$

Define a zero matrix $\mathbb{O}_1=\left[ 0 \right] _{i=1,j=1}^{N_{bp},N_{bp}}$ whose size is $N_{bp}\times N_{bp}$ , $\mathbb{O}_2=\left[ 0 \right] _{i=1,j=1}^{N_{b},N_{bp}}$ whose size is $N_{b}\times N_{bp}$ , $\mathbb{O}_3=\left[ 0 \right] _{i=1,j=1}^{N_{b},N_{b}}$ whose size is $N_{b}\times N_{b}$ , then the block mass matrix is 
$$
M=\left[ \begin{matrix}
	M_e&		\mathbb{O} _3&		\mathbb{O} _2\\
	\mathbb{O} _{3}^{T}&		M_e&		\mathbb{O} _2\\
	\mathbb{O} _{2}^{T}&		\mathbb{O} _{2}^{T}&		\mathbb{O} _1\\
\end{matrix} \right] 
$$

We obtain the first order ODE system
$$
M\vec{X}^{\prime}\left( t \right) +A\vec{X}\left( t \right) =\vec{b}\left( t \right) 
$$

## 8.3 Iteration scheme

$$
\bar{A}\vec{X}^{m+1}=\bar{\vec{b}}^{m+1}, m=0,...,M_m-1
$$

where
$$
\begin{aligned}
&\bar{A}=\frac{M}{\Delta t}+\theta A
\\
&\bar{\vec{b}}^{m+1}=\theta \vec{b}\left( t_{m+1} \right) +\left( 1-\theta \right) \vec{b}\left( t_m \right) +\left[ \frac{M}{\Delta t}-\left( 1-\theta \right) A \right] \vec{X}^m
\end{aligned}
$$

## 8.4 Another format of full discretization

We consider the full discretization (without considering the Dirichlet boundary condition, which will be handled later): for $m=0,...,M_m-1$ , find $\boldsymbol{u}_h^{m+1} \in [U_h]^2$ and $p_h^{m+1}$ in $W_h$ such that
$$
\begin{aligned}
&\left( \frac{\boldsymbol{u}_{h}^{m+1}-\boldsymbol{u}_{h}^{m}}{\Delta t},\boldsymbol{v} \right) +\theta a\left( \boldsymbol{u}_{h}^{m+1},\boldsymbol{v}_h \right) +\left( 1-\theta \right) a\left( \boldsymbol{u}_{h}^{m},\boldsymbol{v}_h \right) 
\\
&+\theta b\left( \boldsymbol{v}_h,p_{h}^{m+1} \right) +\left( 1-\theta \right) b\left( \boldsymbol{v}_h,p_{h}^{m} \right) 
\\
&=\theta \left( \boldsymbol{f}\left( t_{m+1} \right) ,\boldsymbol{v}_h \right) +\left( 1-\theta \right) \left( \boldsymbol{f}\left( t_m \right) ,\boldsymbol{v}_h \right) ,
\\
&\theta b\left( \boldsymbol{u}_{h}^{m+1},q_h \right) +\left( 1-\theta \right) b\left( \boldsymbol{u}_{h}^{m},q_h \right) =0
\end{aligned}
$$

## 8.5 Unsteady linear elasticity equation

Consider
$$
\begin{cases}
	\boldsymbol{u}_{tt}-\nabla \cdot \mathbf{\sigma }\left( \boldsymbol{u} \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \partial \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{u}_0, \frac{\partial \boldsymbol{u}}{\partial t}=\boldsymbol{u}_{00}\,\,\mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega\\
\end{cases}
$$

Define the basic mass matrix
$$
M_e=\left[ m_{ij} \right] _{i,j=1}^{N_b}=\left[ \int_{\Omega}{\phi _j\phi _i\,\,\mathrm{d}x\mathrm{d}y} \right] _{i,j=1}^{N_b}
$$

$\mathbb{O}_3=\left[ 0 \right] _{i=1,j=1}^{N_{b},N_{b}}$ whose size is $N_{b}\times N_{b}$ , then the block mass matrix is 
$$
M=\left[ \begin{matrix}
	M_e&		\mathbb{O} _3\\
	\mathbb{O} _{3}^{T}&		M_e\\
\end{matrix} \right] 
$$

Consider the centered finite difference scheme for the system of ODEs:
$$
M\vec{X}''\left( t \right) +A\vec{X}\left( t \right) =\vec{b}\left( t \right) 
$$

Then the centered finite difference scheme is
$$
\begin{aligned}
&M\frac{\vec{X}^{m+1}-2\vec{X}^m+\vec{X}^{m-1}}{\Delta t^2}+A\frac{\vec{X}^{m+1}+2\vec{X}^m+\vec{X}^{m-1}}{4}
\\
=&\vec{b}\left( t_m \right) ,\ m=1,...,M_m
\end{aligned}
$$

Iteration scheme 2:
$$
\bar{A}\vec{X}^{m+1}=\bar{\vec{b}}^{m+1},\ m=1,...,M_m
$$

where
$$
\begin{aligned}
&\bar{A}=\frac{M}{\Delta t^2}+\frac{A}{4}
\\
&\bar{\vec{b}}^{m+1}=\vec{b}\left( t_m \right) +\left[ \frac{2M}{\Delta t^2}-\frac{A}{2} \right] \vec{X}^m-\left[ \frac{M}{\Delta t^2}+\frac{A}{4} \right] \vec{X}^{m-1}
\end{aligned}
$$

# Chapter 9: Finite elements for 2D unsteady Navier-Stokes equations

## 9.1 Weak formulation

Consider the 2D unsteady Navier-Stokes equation
$$
\begin{cases}
	\boldsymbol{u}_t+\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega \times \left[ 0,T \right]\\
	\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{g}\,\,\mathrm{on}\ \partial \Omega \times \left[ 0,T \right]\\
	\boldsymbol{u}=\boldsymbol{u}_0, p=p_0, \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega\\
\end{cases}
$$

where $\Omega$ is a 2D domain, $[0,T]$ is the time interval, $\boldsymbol{f}(x,y,t)$ is a given function on $\Omega\times[0,T]$ , $\boldsymbol{g}(x,y,t)$ is a given function on $\partial\Omega\times[0,T]$ , $\boldsymbol{u}_0(x,y)$ and $p_0(x,y)$ are given functions on $\Omega$ at $t=0$ , $\boldsymbol{u}(x,y,t)$ and $p(x,y,t)$ are the unknown functions, and
$$
\boldsymbol{u}\left( x,y,t \right) =\left( u_1,u_2 \right) ^T,\boldsymbol{g}\left( x,y,t \right) =\left( g_1,g_2 \right) ^T,\\
\boldsymbol{f}\left( x,y,t \right) =\left( f_1,f_2 \right) ^T,
\boldsymbol{u}_0\left( x,y \right) =\left( u_{10},u_{20} \right) ^T,
$$

Weak formulation
$$
\begin{aligned}
	&\int_{\Omega}{\boldsymbol{u}_t\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{2\nu \mathbb{D} \left( \boldsymbol{u} \right) :\mathbb{D} \left( \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}\\
	&-\int_{\Omega}{p\left( \nabla \cdot \boldsymbol{v} \right) \,\,\mathrm{d}x\mathrm{d}y}-\int_{\partial \Omega}{\left( \mathbb{T} \left( \boldsymbol{u},p \right) \boldsymbol{n} \right) \cdot \boldsymbol{v}\,\,\mathrm{d}s}=\int_{\Omega}{\boldsymbol{f}\cdot \boldsymbol{v}\,\,\mathrm{d}x\mathrm{d}y,}\\
	&-\int_{\Omega}{\left( \nabla \cdot \boldsymbol{u} \right) q\,\,\mathrm{d}x\mathrm{d}y}=0\\
\end{aligned}
$$

## 9.2 Discretization formulation

Semi-discretization
$$
\begin{aligned}
\left( \boldsymbol{u}_{h_t},\boldsymbol{v}_h \right) +c\left( \boldsymbol{u}_h,\boldsymbol{u}_h,\boldsymbol{v}_h \right) +a\left( \boldsymbol{u}_h,\boldsymbol{v}_h \right) +b\left( \boldsymbol{v}_h,p_h \right) &=\left( \boldsymbol{f},\boldsymbol{v}_h \right) 
\\
b\left( \boldsymbol{u}_h,q_h \right) &=0
\end{aligned}
$$

Full discretization
$$
\left( \frac{\boldsymbol{u}_{h}^{m+1,\left( l \right)}-\boldsymbol{u}_{h}^{m}}{\Delta t},\boldsymbol{v}_h \right) 
\\
+\theta \left[ c\left( \boldsymbol{u}_{h}^{m+1,\left( l \right)},\boldsymbol{u}_{h}^{m+1,\left( l-1 \right)},\boldsymbol{v}_h \right) +c\left( \boldsymbol{u}_{h}^{m+1,\left( l-1 \right)},\boldsymbol{u}_{h}^{m+1,\left( l \right)},\boldsymbol{v}_h \right) \right] 
\\
-\theta c\left( \boldsymbol{u}_{h}^{m+1,\left( l-1 \right)},\boldsymbol{u}_{h}^{m+1,\left( l-1 \right)},\boldsymbol{v}_h \right) +\left( 1-\theta \right) c\left( \boldsymbol{u}_{h}^{m},\boldsymbol{u}_{h}^{m},\boldsymbol{v}_h \right) 
\\
+\theta a\left( \boldsymbol{u}_{h}^{m+1,\left( l \right)},\boldsymbol{v}_h \right) +\left( 1-\theta \right) a\left( \boldsymbol{u}_{h}^{m},\boldsymbol{v}_h \right) 
\\
+\theta b\left( \boldsymbol{v}_h,p_{h}^{m+1,\left( l \right)} \right) +\left( 1-\theta \right) \theta b\left( \boldsymbol{v}_h,p_{h}^{m} \right) 
\\
=\theta \left( \boldsymbol{f}\left( t_{m+1} \right) ,\boldsymbol{v}_h \right) +\left( 1-\theta \right) \left( \boldsymbol{f}\left( t_m \right) ,\boldsymbol{v}_h \right) 
\\
\theta b\left( \boldsymbol{u}_{h}^{m+1,\left( l \right)},q_h \right) +\left( 1-\theta \right) b\left( \boldsymbol{u}_{h}^{m},q_h \right) =0
$$

Full discretization in scalar formulation
$$
\int_{\Omega}{\frac{u_{1h}^{m+1,\left( l \right)}-u_{1h}^{m}}{\Delta t}v_{1h}\mathrm{d}x\mathrm{d}y}+\int_{\Omega}{\frac{u_{2h}^{m+1,\left( l \right)}-u_{2h}^{m}}{\Delta t}v_{2h}\mathrm{d}x\mathrm{d}y}
\\
+\theta \int_{\Omega}{\left( u_{1h}^{m+1,\left( l \right)}\frac{\partial u_{1h}^{m+1,\left( l-1 \right)}}{\partial x}v_{1h}+u_{2h}^{m+1,\left( l \right)}\frac{\partial u_{1h}^{m+1,\left( l-1 \right)}}{\partial y}v_{1h}+u_{1h}^{m+1,\left( l \right)}\frac{\partial u_{2h}^{m+1,\left( l-1 \right)}}{\partial x}v_{2h}+u_{2h}^{m+1,\left( l \right)}\frac{\partial u_{2h}^{m+1,\left( l-1 \right)}}{\partial y}v_{2h} \right) \,\,\mathrm{d}x\mathrm{d}y}
\\
+\theta \int_{\Omega}{\left( u_{1h}^{m+1,\left( l-1 \right)}\frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial x}v_{1h}+u_{2h}^{m+1,\left( l-1 \right)}\frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial y}v_{1h}+u_{1h}^{m+1,\left( l-1 \right)}\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial x}v_{2h}+u_{2h}^{m+1,\left( l-1 \right)}\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial y}v_{2h} \right) \,\,\mathrm{d}x\mathrm{d}y}
\\
-\theta \int_{\Omega}{\left( u_{1h}^{m+1,\left( l-1 \right)}\frac{\partial u_{1h}^{m+1,\left( l-1 \right)}}{\partial x}v_{1h}+u_{2h}^{m+1,\left( l-1 \right)}\frac{\partial u_{1h}^{m+1,\left( l-1 \right)}}{\partial y}v_{1h}+u_{1h}^{m+1,\left( l-1 \right)}\frac{\partial u_{2h}^{m+1,\left( l-1 \right)}}{\partial x}v_{2h}+u_{2h}^{m+1,\left( l-1 \right)}\frac{\partial u_{2h}^{m+1,\left( l-1 \right)}}{\partial y}v_{2h} \right) \,\,\mathrm{d}x\mathrm{d}y}
\\
+\left( 1-\theta \right) \int_{\Omega}{\left( u_{1h}^{m}\frac{\partial u_{1h}^{m}}{\partial x}v_{1h}+u_{2h}^{m}\frac{\partial u_{1h}^{m}}{\partial y}v_{1h}+u_{1h}^{m}\frac{\partial u_{2h}^{m}}{\partial x}v_{2h}+u_{2h}^{m}\frac{\partial u_{2h}^{m}}{\partial y}v_{2h} \right) \,\,\mathrm{d}x\mathrm{d}y}
\\
+\theta \int_{\Omega}{\nu \left( 2\frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial x}\frac{\partial v_{1h}}{\partial x}+2\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial y}\frac{\partial v_{2h}}{\partial y}+\frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial y}\frac{\partial v_{1h}}{\partial y}+\frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial y}\frac{\partial v_{2h}}{\partial x}+\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial x}\frac{\partial v_{1h}}{\partial y}+\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial x}\frac{\partial v_{2h}}{\partial x} \right) \,\mathrm{d}x\mathrm{d}y}
\\
+\left( 1-\theta \right) \int_{\Omega}{\nu \left( 2\frac{\partial u_{1h}^{m}}{\partial x}\frac{\partial v_{1h}}{\partial x}+2\frac{\partial u_{2h}^{m}}{\partial y}\frac{\partial v_{2h}}{\partial y}+\frac{\partial u_{1h}^{m}}{\partial y}\frac{\partial v_{1h}}{\partial y}+\frac{\partial u_{1h}^{m}}{\partial y}\frac{\partial v_{2h}}{\partial x}+\frac{\partial u_{2h}^{m}}{\partial x}\frac{\partial v_{1h}}{\partial y}+\frac{\partial u_{2h}^{m}}{\partial x}\frac{\partial v_{2h}}{\partial x} \right) \,\mathrm{d}x\mathrm{d}y}
\\
-\theta \int_{\Omega}{\left( \frac{\partial v_{1h}}{\partial x}p_{h}^{m+1,\left( l \right)}+\frac{\partial v_{2h}}{\partial y}p_{h}^{m+1,\left( l \right)} \right) \mathrm{d}x\mathrm{d}y}-\left( 1-\theta \right) \int_{\Omega}{\left( \frac{\partial v_{1h}}{\partial x}p_{h}^{m}+\frac{\partial v_{2h}}{\partial y}p_{h}^{m} \right) \mathrm{d}x\mathrm{d}y}
\\
=\theta \int_{\Omega}{\left( f_1\left( t_{m+1} \right) v_{1h}+f_2\left( t_{m+1} \right) v_{2h} \right) \mathrm{d}x\mathrm{d}y}+\left( 1-\theta \right) \int_{\Omega}{\left( f_1\left( t_m \right) v_{1h}+f_2\left( t_m \right) v_{2h} \right) \mathrm{d}x\mathrm{d}y}
\\
-\theta \int_{\Omega}{\left( \frac{\partial u_{1h}^{m+1,\left( l \right)}}{\partial x}q_h+\frac{\partial u_{2h}^{m+1,\left( l \right)}}{\partial y}q_h \right) \mathrm{d}x\mathrm{d}y}-\left( 1-\theta \right) \int_{\Omega}{\left( \frac{\partial u_{1h}^{m}}{\partial x}q_h+\frac{\partial u_{2h}^{m}}{\partial y}q_h \right) \mathrm{d}x\mathrm{d}y}=0
$$

## 9.3 Matrix formulation

$$
\left( \frac{M}{\Delta t}+\theta A+\theta AN^{m+1} \right) X^{m+1,\left( l \right)}=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) +\theta \overrightarrow{bN^{m+1}}-\left( 1-\theta \right) \overrightarrow{bN^m}+\left[ \frac{M}{\Delta t}-\left( 1-\theta \right) A \right] X^m
$$

Main pseudo code:
- Generate the mesh information matrices P and T
- Assemble the mass matrix M and stiffness matrix A
- Generate the initial vector $\vec{X}^0$
- Iterate in time: $FOR\ m=0,...,M_m-1$
  - $t_{m+1}=(m+1) \Delta t$
  - Assemble the load vector $b(t_{m+1})$
  - Newton iteration: $FOR\ l=1,2,...,L$
    - Assemble the matrix $AN$
    - Assemble the vector $bN$
    - $A^{m+1,(l)}$ and $\vec{b}^{m+1,(l)}$
    - Treat Dirichlet boundary condition
    - Solve $A^{m+1,(l)}\vec{X}^{m+1,(l)}=\vec{b}^{m+1,(l)}$
  - $END$
  - Let $\vec{X}^{m+1}$ be the final $\vec{X}^{m+1,(l)}$ from the Newton's iteration
- $END$

## 9.4 Numerical example

Example 1: On the domain $\Omega = [0, 1]\times[−0.25, 0]$, consider
the time-dependent Navier-Stokes equation
$$
\begin{aligned}
&\boldsymbol{u}_t+\left( \boldsymbol{u}\cdot \nabla \right) \boldsymbol{u}-\nabla \cdot \mathbb{T} \left( \boldsymbol{u},p \right) =\boldsymbol{f}\,\,\mathrm{in}\ \Omega \times \left[ 0,T \right] 
\\
&\nabla \cdot \boldsymbol{u}=0\ \mathrm{in}\ \Omega \times \left[ 0,T \right] 
\\
&u_1=x^2y^2+e^{-y},\ \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega 
\\
&u_2=-\frac{2}{3}xy^3+2-\pi \sin \left( \pi x \right) ,\ \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega 
\\
&p=-\left[ 2-\pi \sin \left( \pi x \right) \right] \cos \left( 2\pi y \right),\ \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega
\end{aligned}
$$

Continued formulation:
$$
\begin{aligned}
&u_1=e^{-y}\cos \left( 2\pi t \right) \,\,\mathrm{on}\ x=0
\\
&u_1=\left( y^2+e^{-y} \right) \cos \left( 2\pi t \right) \,\,\mathrm{on}\ x=1
\\
&u_1=\left( \frac{1}{16}x^2+e^{0.25} \right) \cos \left( 2\pi t \right) \,\,\mathrm{on}\ y=-0.25
\\
&u_1=\cos \left( 2\pi t \right) \,\,\mathrm{on}\ y=0
\\
&u_2=2\cos \left( 2\pi t \right) \,\,\mathrm{on}\ x=0
\\
&u_2=\left( -\frac{2}{3}y^3+2 \right) \cos \left( 2\pi t \right) \,\,\mathrm{on}\ x=1
\\
&u_2=\left[ \frac{1}{96}x+2-\pi \sin \left( \pi x \right) \right] \cos \left( 2\pi t \right) \,\,\mathrm{on}\ y=-0.25
\\
&u_2=\left[ 2-\pi \sin \left( \pi x \right) \right] \cos \left( 2\pi t \right) \,\,\mathrm{on}\ y=0
\end{aligned}
$$

Here
$$
\begin{aligned}
f_1=&-2\pi \left( x^2y^2+e^{-y} \right) \sin \left( 2\pi t \right) 
\\
&+\left\{ 2xy^2\left( x^2y^2+e^{-y} \right) +\left[ -\frac{2}{3}xy^3+2-\pi \sin \left( \pi x \right) \right] \left( 2yx^2-e^{-y} \right) \right\} \cos ^2\left( 2\pi t \right) 
\\
&+\left[ -2\nu x^2-2\nu y^2-\nu e^{-y}+\pi ^2\cos \left( \pi x \right) \cos \left( 2\pi y \right) \right] \cos \left( 2\pi t \right) 
\\
f_2=&-2\pi \left[ -\frac{2}{3}xy^3+2-\pi \sin \left( \pi x \right) \right] \sin \left( 2\pi t \right) 
\\
&+\left\{ \left( x^2y^2+e^{-y} \right) \left[ -\frac{2}{3}y^3-\pi ^2\cos \left( \pi x \right) \right] -2xy^2\left[ -\frac{2}{3}xy^3+2-\pi \sin \left( \pi x \right) \right] \right\} \cos ^2\left( 2\pi t \right) 
\\
&+\left[ 4\nu xy-\nu \pi ^3\sin \left( \pi x \right) +2\pi \left( 2-\pi \sin \left( \pi x \right) \right) \sin \left( 2\pi y \right) \right] \cos \left( 2\pi t \right) 
\end{aligned}
$$

The analytic solution of this problem is
$$
\begin{aligned}
&u_1=\left( x^2y^2+e^{-y} \right) \cos \left( 2\pi t \right) 
\\
&u_2=\left[ -\frac{2}{3}xy^3+2-\pi \sin \left( \pi x \right) \right] \cos \left( 2\pi t \right) 
\\
&p=-\left[ 2-\pi \sin \left( \pi x \right) \right] \cos \left( 2\pi y \right) \cos \left( 2\pi t \right) 
\end{aligned}
$$

# 进一步学习

## 程序优化

- Dirichlet 边界条件处理所用时间太长，是否有其他计算方法(可以将目前这种计算方法保留，添加一个新的输入参数表示使用哪种计算方法即可)
- 优化当前计算所使用的脚本，将前处理、求解器、后处理进一步打包分离
- 将与课件可以进行对比的算例制作成标准测试算例
- 目前的做法是，将不同方程的组装方法写成了多个类，但实际上这些类有许多通用的地方，是否可以将这些类整合起来，进一步优化程序结构
- 求解类型也可以想办法整合，因为有稳态求解和非稳态求解，看看怎么合并
- 方程组装中涉及到了稀疏矩阵的分块组装，自己在里面写了很多重复的操作，并且每一次写不同的程序还需要思考组装的索引等这些琐碎的细节，是否可以写一个同的分块组装子程序
- 管理程序文件，分好类，如 math 文件夹、element 文件夹、equation 文件夹等

## Slides 拓展内容整理

- 第二章
  - 四边形单元
  - 3D单元
- 第三章
  - mixed boundary conditions
  - Non-isotropic equation
  - A more general second order equation
- 第四章
  - Mixed boundary conditions
  - Non-isotropic second order parabolic equation with mixed
boundary conditions
  - Second order hyperbolic equation
  - Non-isotropic second order hyperbolic equation with mixed
boundary conditions
- 第五章
  - stress, robin, dirichlet/stress/robin mixed boundary conditions
  - stress, robin, dirichlet/stress/robin mixed boundary conditions in normal/tangential directions
- 第六章
  - stress, robin, dirichlet/stress/robin mixed boundary conditions
  - stress, robin, dirichlet/stress/robin mixed boundary conditions in normal/tangential directions
- 第七章
  - stress, robin, dirichlet/stress/robin mixed boundary conditions
  - stress, robin, dirichlet/stress/robin mixed boundary conditions in normal/tangential directions
- 第八章
  - mixed boundary conditions
  - mixed boundary conditions in normal/tangential directions
    - 非稳态条件下的 stress, robin boundary conditions 是否需要在每一个时间迭代步都进行边界条件处理？
    - 第四章视频下的评论：非稳态问题 $\theta$ 格式的Neumman和Robin边界需要对 $A(t_{m+1})$、$A(t_m)$、$b(t_{m+1})$、$b(t_m)$ 分别处理一遍，否则精度会降低一阶
    - 第四章拓展部分问题，写出完整的弱形式后再进行有限元离散和时间离散，看看到底是不是需要
    - 根据课程视频，如果边界条件与时间相关，需要对 $A(t_{m+1}), A(t_m), b(t_{m+1}), b(t_m)$ 分别处理
    - If the functions in the stress/Robin boundary conditions depend on time, then the same algorithms as those in Chapter 7 can be used at each time iteration step. But the time needs to be specified in these algorithms.
  - unsteady linear elasticity equation
    - 给出初始时刻的边界条件一阶时间偏导，如何给出初始的 $X_1$ ？
    - 目前可能的做法是使用中心差分格式 $2\Delta t X_{0}^{\prime} = X_1 - X_{-1}$ 和递推格式 $AX_1=b_1$ 一起求解出 $X_1$
    - mixed boundary conditions
    - mixed boundary conditions in normal/tangential directions
- 第九章
  - mixed boundary conditions
  - mixed boundary conditions in normal/tangential directions