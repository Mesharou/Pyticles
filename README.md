Version of Python 3.6 running

options: Tested

x_periodic = False
y_periodic = False
adv3d = True
meanflow=False
timing = True
write_uv=False
restart = False
initial_depth = False
continuous_injection = False

======================================
Python portability to python 3 In bugs : 

- numpy retruns a float for integer division,
  the solution is to use the operator //
ex :
10 / 2 
> 5.0 

10 // 2
> 5
---------------------------------------





