## Dec 9, 2018

添加可成功运行的CUDA代码，结果与c代码相同，性能待优化

### 编译:

加`-run`选项在编译后立即运行

```sh
nvcc -o FDTD2Dsolver FDTD2Dsolver.cu [-run]
```

### 运行、性能优化

```sh
nvprof ./FDTD2Dsolver
```

输出示例：

```
==1270== Profiling application: ./FDTD2Dsolver
==1270== Profiling result:
            Type  Time(%)      Time     Calls       Avg       Min       Max  Name
 GPU activities:   76.25%  35.085ms       500  70.169us  67.999us  73.855us  iteration(float*, float*, float, float, int, int)
                   23.45%  10.791ms         2  5.3954ms  5.3673ms  5.4235ms  [CUDA memcpy DtoH]
                    0.18%  82.015us         1  82.015us  82.015us  82.015us  exactSolution(float*, int, int, float, float, float)
                    0.12%  55.231us         1  55.231us  55.231us  55.231us  initArray(float*, int, int, float, float)
      API calls:   71.11%  139.10ms         2  69.549ms  183.45us  138.92ms  cudaMalloc
                   19.03%  37.231ms       501  74.313us  59.347us  93.410us  cudaDeviceSynchronize
                    5.69%  11.122ms         2  5.5611ms  5.5280ms  5.5943ms  cudaMemcpy
                    2.26%  4.4225ms       502  8.8090us  8.0920us  31.747us  cudaLaunch
                    1.02%  1.9867ms      3011     659ns     304ns  14.513us  cudaSetupArgument
                    0.44%  858.88us         2  429.44us  229.48us  629.40us  cudaFree
                    0.18%  353.01us       502     703ns     659ns  2.9670us  cudaConfigureCall
                    0.14%  270.68us        94  2.8790us     759ns  70.136us  cuDeviceGetAttribute
                    0.13%  251.30us         1  251.30us  251.30us  251.30us  cuDeviceTotalMem
                    0.01%  17.952us         1  17.952us  17.952us  17.952us  cuDeviceGetName
                    0.00%  3.6320us         3  1.2100us     756ns  1.8470us  cuDeviceGetCount
                    0.00%  2.0340us         2  1.0170us     818ns  1.2160us  cuDeviceGet
```

生成图像同下

```sh
python plot.py
```

---

## Nov 27, 2018

当前阶段：直接从MATLAB代码转化为c

### 编译:

```sh
gcc -Wall -g -o FDTD2Dsolver FDTD2Dsolver.c -lm
```

### 运行:

```sh
./FDTD2Dsolver
python plot.py
```

### 结果：

生成3个数据文件`u0.txt`, `exact.txt`, `diff.txt`

python脚本利用上述数据生成图像


vx=vy=0.5,NX=NY=2048,Ntimesteps=500,结束时间约t=19.5

![Result](plot.png)
