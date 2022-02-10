# sage code in C++
## BUPT RTT2 Lab
Sage algorithm for Channel Modeling


## 简易流程

#### 1. 准备数据

- 首先在/home/rtt2/code/sage/data下新建自己的文件夹：xxx
- 将自己的的CIR bin文件以spot为文件夹放到目录：/home/rtt2/code/sage/data/xxx/dd,此处dd目前应该必须是99以内数字
- 将自己测量使用的天线校准文件放到目录：/home/rtt2/code/sage/data/xxx/antenna_bin
- 同样在同一级目录下新建result文件夹：/home/rtt2/code/sage/data/xxx/result

*****

#### 2. 修改代码

**主要需要修改的代码只有两个文件：test_mpi_sage.cpp(这个是运行主文件)，isi_sage.h(主要是修改参数配置的地方，后续直接暴露接口传参，不用改代码)**

1. 在isi_sage.h 修改ConfigInit的参数列表中的几个参数
2. 在test_mpi_sage.cpp的面函数里面,主要修改spots , 以及对应的三个前面要求建立的路径，另外还要修改文件开头的ParseFileNum函数中cycle的起始序号和终止序号

*****
#### 3. 编译代码

**分布式并行编译代码需要按照以下步骤在make文件同级目录下执行：**

```bash
# 1. 清理掉以前编译的文件
make clean

# 2. 重新编译
make

# 3. 运行程序
nohup mpiexec -f hostfile -n 50 ./test_mpi_sage > sage_20181220.log 2>&1 &
## 解释以上代码：
## nohup: linux 后台静默执行程序指令
## mpiexec: MPICH框架下的程序执行指令
## hostfile: 指定哪几台linux主机作为分布式计算机器，在make同级目录下，
## hostfile中的机器需要提前设置对应机器的ssh免密访问，后续会有教程介绍部署过程，目前单机下需要使用注释掉后面几台机器，只保留第一行主机信息
## -n 表示并行度，必须要小于所有机器的总CPU核数，建议每台机器使用CPU核数的80%,否则程序一跑，机器就会无响应，然后挂掉
## > sage_20181220.log 表示将打印的信息重定向到log文件，自己取名字，系统会自动建立log文件
## 2>&1 & 表示将错误信息重定向标准输出流，标准输出流已经重定向到了log文件，因此错误信息也会打印到log文件
## 必须要打印到log文件，一个是调试，另外一个是分布式计算是并行的，log信息比较多，直接运行命令行终端会一直打印信息，没办法干其他事

```
如果前面都正常执行，最后程序执行完毕之后会在result文件下产生以cycle为单位的结果文件，为二进制文件，后续可以改为二进制和文本文件可选


### 后续会完善教程。。。





