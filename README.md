# Домашнее задание: Прямой и Дуальный Симплекс Метод

Я реализовал оба метода в рамках одного файла. В def Solve(c, A, b, mode) можно подать mode=primal для обычного, dual для дуального или auto для автоопределения.

```b
<details><summary>Вывод программы</summary><p>

<pre>
$ python simplex_template <example_phase1.txt>
Running SIMPLEX.py with file .\example_phase1.txt...

n = 2
m = 3
c = [2. 1.]
A = [[-1.  1.]
 [-1. -2.]
 [ 0.  1.]]
b = [-1. -2.  1.]
Solving the linear program using the Primal Simplex method...


Result:
Status: unbounded
The problem is unbounded.

$ python simplex_template <example_phase1.txt>
Running SIMPLEX.py with file .\example_stalling.txt...

n = 4
m = 3
c = [ 1. -2.  0. -2.]
A = [[ 0.5 -3.5 -2.   4. ]
 [ 0.5  0.  -0.5  0.5]
 [ 1.   0.   0.   0. ]]
b = [0. 0. 1.]
Solving the linear program using the Primal Simplex method...


Result:
Status: optimal
Optimal solution x* = [1. 0. 1. 0.]
Optimal value = 1.0
</pre>

</p></details>
```
