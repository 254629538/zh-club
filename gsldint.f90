Module gsld !高斯—勒让德积分高斯点及权重的求解模块
  Implicit None
  Integer,Parameter :: nk = 8                          !设置求解高斯点的个数
  !Integer,Parameter :: DP = Selected_Real_Kind( p=13 ) !设置kind的数值
  Real*8, parameter :: EPS = 1.0D-15          !精度设置
  real*8::fn1(nk),ak1(nk)
  save fn1,ak1

Contains

  Real*8 Function f(x) !定义函数f（x）
    Implicit None
    Integer :: i
    Real*8 :: a(nk), x !a(n)代表n阶勒让德多项式
    a(1) = x !1阶勒让德多项式
    a(2) = 1.5D0*(x**2) - 0.5D0 !2阶勒让德多项式
    Do i = 3, nk
      a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
      !利用递推关系产生n阶勒让德多项式
    End Do
    f = a(nk) !生成的n阶勒让德多项式
  End Function f

  Real*8 Function f1(x)
    Implicit None
    Integer :: i
    Real*8 :: a(nk), x
    a(1) = x
    a(2) = 1.5D0*x**2 - 0.5D0
    Do i = 3, nk - 1
      a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
    End Do
    f1 = a(nk-1) !生成的（n-1）阶勒让德多项式
  End Function f1

  Real*8 Function g(x)
    Implicit None
    Integer :: i
    Real*8 :: a(nk), x
    a(1) = x
    a(2) = 1.5D0*x**2 - 0.5D0
    Do i = 3, nk
      a(i) = (2*i-1)*x*a(i-1)/i - (i-1)*a(i-2)/i
    End Do
    g = nk*a(nk-1)/(1-x**2) - nk*x*a(nk)/(1-x**2)
    !生成n阶勒让德多项式的导数表达式
  End Function g

  Real*8 Function bis(a, b) !二分法求解函数的解
    Implicit None
    Real*8:: a, b, c
    !a,b是传递进来的划分好的有一个解存在的区间
    Do
      c = (a+b)/2.0D0
      If (f(c)*f(a)<0) Then
        b = c
      Else
        a = c
      End If
      If ((b-a)<EPS) exit
    End Do
    bis = c!bis即是利用二分法求得的解
  End Function bis

  Subroutine fn0(fn, ak)
    Implicit None
    Real*8 :: m, fn(nk), ak(nk)
    !定义数组,大小n由module开始声明。
    Integer :: i, j
    j = 0 !赋值控制循环变量的初值
    m = -1.000001 !设置计算域[-1，1] 的下限，即代替-1
    Do i = 1, 200000 !这个循环次数应该是由步长0.00001决 定,计算方法：200000=2/0.000001
      If (f(m)*f(m+0.00001)<0) Then !从下限处开始往上逐步累加，
        !由步长0.00001说明最多求解10^5个解
        j = j + 1 !记录这是第几个解
        fn(j) = bis(m, m+0.00001)
        !调用二分法求解程序在分好的一小段上求解，
        !将解存储在fn（j）
        ak(j) = 2.0D0/(nk*f1(fn(j))*g(fn(j))) !高斯点的权重
        ! Write (*, *) '高斯点序号', j
        ! Write (*, *) '高斯点', fn(j)
        ! Write (*, *) '高斯点权重', ak(j)
      End If
      m = m + 0.00001 !执行完一次判断m向前推进一步
    End Do
  End Subroutine fn0

End Module gsld
