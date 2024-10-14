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
\\
u\left( a \right) =g_a, \  u\left( b \right) =g_b
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
\\
E_n=\left[ x_n,x_{n+1} \right] \,\,\left( n=1,\cdots ,N \right) 
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
\\
u_h=\sum_{j=1}^{N+1}{u_j\phi _j}
\end{gather}
$$
将test function选为同样的分片线性函数 $v_h=\phi_i \  \left( i=1,\cdots,N+1 \right)$ 
> 因为 $u_j \ (j=1,\cdots,N+1)$ 共 $(N+1)$ 个未知数需要求解，所以选择 $(N+1)$ 个test function测试，最后形成方程组联立求解


$$
\begin{align}
&\int_a^b{c\left( \sum_{j=1}^{N+1}{u_j\phi _j} \right) ^{\prime}\phi _{i}^{\prime}\mathrm{d}x=\int_a^b{f\phi _{i}\mathrm{d}x}}
,\quad i=1,\cdots,N+1
\\
\Rightarrow &\sum_{j=1}^{N+1}{\left[ \int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x} \right]u_j}=\int_a^b{f\phi _{i}\mathrm{d}x}
,\quad i=1,\cdots,N+1
\end{align}
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
\\
b_i=\int_a^b{f\phi _{i}\mathrm{d}x}
\end{gather}
$$
### Assembly
$$
\begin{gather}
A_{ij}=\int_a^b{c\phi _{j}^{\prime}\phi _{i}^{\prime}\mathrm{d}x}=\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{c\psi _{nj}^{\prime}\psi _{ni}^{\prime}\mathrm{d}x}},i,j=1,\cdots ,N+1
\\
b_i=\int_a^b{f\phi _i\mathrm{d}x}=\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{f\psi _{ni}\mathrm{d}x}},i,j=1,\cdots ,N+1
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
\\
\frac{\mathrm{d}\psi _{nj}\left( x \right)}{\mathrm{d}x}=\frac{\mathrm{d}\hat{\psi}_j\left( \hat{x} \right)}{\mathrm{d}\hat{x}}\frac{\mathrm{d}\hat{x}}{\mathrm{d}x}
\end{gather}
$$

### Neumann/Robin boundary conditions

#### Neumann
$$
\begin{gather}
-\frac{\mathrm{d}}{\mathrm{d}x}\left( c\left( x \right) \frac{\mathrm{d}u\left( x \right)}{\mathrm{d}x} \right) =f\left( x \right) , a<x<b
\\
u^{\prime}\left( a \right) =r_a, u\left( b \right) =g_b
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
\\
u\left( a \right) =g_a,u^{\prime}\left( b \right) +q_bu\left( b \right) =p_b
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
\begin{align}
\left\| u-u_h \right\| _{\infty}&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -u_h\left( x \right) \right|
\\
&=\,\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -\sum_{j=1}^{N_b}{u_j\phi _j\left( x \right)} \right|
\\
&=\underset{1\le n\le N}{\max}\underset{x_n\le x\le x_{n+1}}{\max}\left| u\left( x \right) -\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}\left( x \right)} \right|
\end{align}
$$
$$
\begin{align}
\left\| u-u_h \right\| _0&=\sqrt{\int_I{\left( u-u_h \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-u_h \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{j=1}^{N_b}{u_j\phi _j} \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}} \right) ^2\mathrm{d}x}}}
\end{align}
$$
$$
\begin{align}
\left| u-u_h \right|_1&=\sqrt{\int_I{\left( u^{\prime}-u_{h}^{\prime} \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u^{\prime}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}^{\prime}} \right) ^2\mathrm{d}x}}}
\end{align}
$$
可以发现，误差计算的关键在于计算下列公式的子程序
$$
w_{n,s}=\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}^{\left( s \right)}}
$$
由于基函数库已经写好了，就只需要理解如何从已经计算好的有限元解 $u_h$ 中取出单元上的有限元局部解向量 $u_{T_b\left( k,n \right)}$ ，这里就不剧透了。

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
\begin{align}
\hat{x}&=\frac{\left( y_3-y_1 \right) \left( x-x_1 \right) -\left( x_3-x_1 \right) \left( y-y_1 \right)}{\left| J \right|}
\\
\hat{y}&=\frac{-\left( y_2-y_1 \right) \left( x-x_1 \right) +\left( x_2-x_1 \right) \left( y-y_1 \right)}{\left| J \right|}
\end{align}
$$
这样可以使用参考基函数
线性基函数
$$
\begin{align}
\hat{\psi}_1\left( \hat{x},\hat{y} \right) &=-\hat{x}-\hat{y}+1
\\
\hat{\psi}_2\left( \hat{x},\hat{y} \right) &=\hat{x}
\\
\hat{\psi}_3\left( \hat{x},\hat{y} \right) &=\hat{y}
\end{align}
$$
二阶基函数
$$
\begin{align}
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
\end{align}
$$
采用链式求导法则
$$
\begin{align}
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
\end{align}
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
\begin{align}
&-\nabla \cdot \left( c\nabla u \right) =f, \ \mathrm{in} \ \Omega 
\\
&u=g,\ \mathrm{on} \  \partial \Omega 
\end{align}
$$
其中梯度和散度为
$$
\begin{align}
&\nabla u=\left( u_x,u_y \right) 
\\
&\nabla \cdot \vec{v}=\frac{\partial v_1}{\partial x}+\frac{\partial v_2}{\partial y}
\end{align}
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
\begin{align}
A_{ij}&=\int_{\Omega}{c\nabla \phi _j\cdot \nabla \phi _i\mathrm{d}x\mathrm{d}y}
\\
b_i&=\int_{\Omega}{f\phi _i\mathrm{d}x\mathrm{d}y}
\end{align}
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
\begin{align}
\left\| u-u_h \right\| _0&=\sqrt{\int_I{\left( u-u_h \right) ^2\mathrm{d}x}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-u_h \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{j=1}^{N_b}{u_j\phi _j} \right) ^2\mathrm{d}x}}}
\\
&=\sqrt{\sum_{n=1}^N{\int_{x_n}^{x_{n+1}}{\left( u-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\psi _{nk}} \right) ^2\mathrm{d}x}}}
\end{align}
$$
$$
\begin{align}
\left| u-u_h \right|_{1}^{2}=&\sum_{n=1}^N{\int_{E_n}{\left( \frac{\partial u}{\partial x}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\frac{\partial \psi _{nk}}{\partial x}} \right) ^2\mathrm{d}x\mathrm{d}y}}
\\
&+\sum_{n=1}^N{\int_{E_n}{\left( \frac{\partial u}{\partial y}-\sum_{k=1}^{N_{lb}}{u_{T_b\left( k,n \right)}\frac{\partial \psi _{nk}}{\partial y}} \right) ^2\mathrm{d}x\mathrm{d}y}}
\end{align}
$$

## 3.5 More Discussion

### 3.5.1 Neumann boundary conditions

$$
\begin{align}
\int_{\partial \Omega}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}&=\int_{\Gamma _N}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}+\int_{\partial \Omega /\Gamma _N}{\left( c\nabla u\cdot \vec{n} \right) v\mathrm{d}s}
\\
&=\int_{\Gamma _N}{cpv\mathrm{d}s}
\end{align}
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

这里需要强调一点，边界边信息和边界节点信息分属于不同类别，边界边为网格概念，边界节点为有限元概念。我在上一章的作业中混淆了这两个概念，刚好后面何老师在答疑课中也讲了这个问题，具体链接如下，时间大约在36:39
[有限元基础编程课程答疑-何晓明-2021-1-07](https://www.bilibili.com/video/BV1ET4y1N77h/?share_source=copy_web&vd_source=4846b273dd993443973d15943a9b6546)

# Chapter 4: Finite Elements for 2D second order parabolic and hyperbolic equation

## 4.1 Weak formulation /4.2 Semi-discretization

Cosider the 2D second order parabolic equation
$$
\begin{align}
&u_t-\nabla \cdot \left( c\nabla u \right) =f, \  \mathrm{in} \ \Omega \times \left[ 0,T \right] 
\\
&u=g,\ \mathrm{on}\ \partial \Omega \times \left[ 0,T \right] 
\\
&u=u_0,\ \mathrm{at}\ t=0\ \mathrm{and}\ \mathrm{in}\ \Omega 
\end{align}
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
\begin{align}
\tilde{A}_{m+1}&=\frac{M}{\Delta t}+\theta A\left( t_{m+1} \right) 
\\
\tilde{b}_{m+1}&=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) +\left( \frac{M}{\Delta t}-\left( 1-\theta \right) A\left( t_m \right) \right) u_m
\end{align}
$$
简单起见，考虑LTI系统
$$
\begin{align}
\tilde{A}&=\frac{M}{\Delta t}+\theta A
\\
\tilde{b}_{m+1}&=\theta b\left( t_{m+1} \right) +\left( 1-\theta \right) b\left( t_m \right) +\left[ \frac{M}{\Delta t}-\left( 1-\theta \right) A \right] u_m
\end{align}
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
\begin{align}
&\int_{\Omega}{ ( \lambda \frac{\partial u_1}{\partial x_1}\frac{\partial v_1}{\partial x_1}+2\mu \frac{\partial u_1}{\partial x_1}\frac{\partial v_1}{\partial x_1}+\lambda \frac{\partial u_2}{\partial x_2}\frac{\partial v_1}{\partial x_1}\,\,}
\\
&+\mu \frac{\partial u_1}{\partial x_2}\frac{\partial v_1}{\partial x_2}+\mu \frac{\partial u_2}{\partial x_1}\frac{\partial v_1}{\partial x_2}+\mu \frac{\partial u_1}{\partial x_2}\frac{\partial v_2}{\partial x_1}+\mu \frac{\partial u_2}{\partial x_1}\frac{\partial v_2}{\partial x_1}
\\
&+\lambda \frac{\partial u_1}{\partial x_1}\frac{\partial v_2}{\partial x_2}+\lambda \frac{\partial u_2}{\partial x_2}\frac{\partial v_2}{\partial x_2}+2\mu \frac{\partial u_2}{\partial x_2}\frac{\partial v_2}{\partial x_2}  )\,\,\mathrm{d}x_1\mathrm{d}x_2
\\
=&\int_{\Omega}{\left( f_1v_1+f_2v_2 \right) \,\, \mathrm{d}x_1\mathrm{d}x_2}
\end{align}
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
\begin{align}
&\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}+2\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&=\int_{\Omega}{f_1\phi _i\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{align}
$$
再取 $\boldsymbol{v}_h=\left( 0,\phi _i \right) ^T$ 
$$
\begin{align}
&\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}+\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{1j}\frac{\partial \phi _j}{\partial x_1}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&+\int_{\Omega}{\lambda \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}+2\int_{\Omega}{\mu \left( \sum_{j=1}^{N_b}{u_{2j}\frac{\partial \phi _j}{\partial x_2}} \right) \frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
&=\int_{\Omega}{f_2\phi _i\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{align}
$$
排列成矩阵形式
$$
\begin{align}
A_1=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2},\quad    A_2=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_3=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2},  \quad   A_4=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_5=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2},    \quad A_6=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_1}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\\
A_7=\int_{\Omega}{\mu \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_1}\,\,\mathrm{d}x_1\mathrm{d}x_2},    \quad A_8=\int_{\Omega}{\lambda \frac{\partial \phi _j}{\partial x_2}\frac{\partial \phi _i}{\partial x_2}\,\,\mathrm{d}x_1\mathrm{d}x_2}
\end{align}
$$
写为，也就是说，只需要多次调用矩阵组装的子程序，最后排成大矩阵求解即可
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

# 学习小结

在最简单的一维问题中打好基础，理解好有限元通用结构后，后续的通用程序大部分部分基本不变，每次只需要做少量的修改。然后逐渐向不同的有限元、不同的边界条件、二维三维、非稳态、非线性情况去拓展自己写出来的程序包。

目前完成的部分算是《有限元基础编程I》，差不多是正常教学一学期的内容，后续的内容《有限元基础编程II》录播中也有，之后打算进一步学习。同时，前半段的学习算是手把手带上路，包括问题描述、公式推导、伪代码课件中都有详细的讲述，指导性编程也比较详细，之后的学习应尝试自己独立走完整个流程，完成课件中提到的一些小课题等内容，自行设计和完成相应的算例测试。

## Review

重新把之前的视频、实现代码、学习记录看一遍，然后总结、重构自己的程序，再开始推进下一阶段的学习。
- ref -> local -> global 仿射变换通用处理
- 边界条件程序化处理，需要更多的信息！！！
- 误差计算程序
- 注意 Slides 中明显的重要模块、计算公式！！！
- 三角形单元和四边形单元实现，线性和二次单元

# Chapter 6: 

