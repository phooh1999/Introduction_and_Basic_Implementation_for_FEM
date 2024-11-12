# 概述
- 天元数学东北中心短期课程：《有限元基础编程》
- 何晓明老师 [MISSOURI S&T Dr.Xiaoming He](https://sites.mst.edu/xiaominghe/)
- [课程录播](https://www.bilibili.com/video/BV1Zv411t7Lj/)
- [优质笔记](https://github.com/chichuDlong/FiniteElementMethod_Book)

**课程主要目标**  
1. 讲解有限元算法通用结构以及相应的基本有限元程序的通用结构；
2. 从零开始逐步写出自己的通用结构有限元程序包（而不是理解和使用某个已有程序包）；
3. 从而掌握给各类方程独立写基本有限元通用结构程序的方法；
4. 具备向功能更强大结构更复杂的有限元程序包扩展的能力。

**仓库内容**
- `Homework` : 课程提供的八次作业、参考程序，以及自己的作业程序
- `MyCode` : 尝试对自己的作业程序进行重构
  - 将局部基函数的索引循环改成向量乘法
  - 计算单元刚度阵，使用三元组组装总刚度阵（稀疏矩阵）
- `Slides` : 课程资料，以及一本不错的[笔记](https://github.com/chichuDlong/FiniteElementMethod_Book)
- `notes` : 课程学习过程中记录的笔记

# Course Q&A
（下面的这些链接好像都访问受限了...）

0. [对这门课的一些进一步解释](https://shimo.im/docs/3grh3gkqD6h39xkh/)
1. [关于课件中格式推导和解释的问题](https://shimo.im/docs/yPgytPVYJkTPKW66/)
2. [关于课件中伪代码，模块化，和参数化的问题](https://shimo.im/docs/39XQycqHwr6G66X8/) 
3. [关于网格剖分及其与有限元核心程序接口的问题](https://shimo.im/docs/rqp9DJ6dr3YyKDYc/)
4. [关于有限元基函数程序模块的问题](https://shimo.im/docs/pJ98vP9KJxdjtTHT/) 
5. [关于矩阵和向量组装程序模块的问题](https://shimo.im/docs/CTHw8jp9DvDGqt9H/) 
6. [关于边界处理程序模块的问题](https://shimo.im/docs/XCpp8xpGCyJvVkGG/) 
7. [关于矩阵计算及其与有限元核心程序接口的问题](https://shimo.im/docs/V8prwr6tt9r9D399/) 
8. [关于算例误差计算的程序模块的问题](https://shimo.im/docs/8kdp3kvjJT9xCYVh/)
9. [关于升级到各类更复杂方程的问题（高维，非稳态，耦合，非线性方程等](https://shimo.im/docs/TQRdhqQdtkTxwYD3/)
10. [其他与本课程相关的问题（例如基础程序学习方法，益处，将来的用法，等等](https://shimo.im/docs/wXxytCxcPW8yTxxR/)