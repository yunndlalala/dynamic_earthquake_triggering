### 安装
1. 将Detect_DynamicTriggering程序包放置于任意路径下
2. 运行setup.py

### 模块介绍
dyntrigger包中共包含四个子模块，分别是betatest，calstress，hifi和utils
#### betatest
计算beta值
#### calstress
计算最大动态应力
#### hifi
HiFi方法检测动态触发
> Step1: 生成功率积分数据库

运行Test/bin/cal_PIdatabase.py脚本，对不同的台站建立不同的数据库
* 输入
    - 台站名称列表 sta_list
    - 要计算的三分量列表 chn_list 
    - 仪器响应函数的相关参数
        * sensivity PZ文件的第22行
        * normalizing PZ文件的第23行
        * zeros 零点
        * poles 极点
    - 计算功率谱时截取的时间段长度 time_segment
    - 积分的频率区间 f_win_list
    - 原始数据的存放路径 data_path 
        * 原始数据的路径必须是按照年/日/某台站的sac文件，例如data_path/2006/20060101/GDXB.BHZ
        * sac文件的命名方式必须是station_name.chn, 例如GDXB.BHZ。参数“台站名称”要与station_name完全相同，且可以包含台网名。
    - 存放结果的路径 out_path
    - 并行计算的进程数 cores
        * 本脚本可进行并行计算，如果cores设置为大于1，则多进程并行计算；如果为1，则单进程计算。

* 输出
    - 每个台站的数据库单独为一个文件夹
    - 文件夹中每天的结果存储为一个csv文件，文件以station_name.chn_day.csv形式命名。例如NC.GDXB.HHZ_20060618.csv
    - csv文件中每行记录一个时间段内不同频率范围的积分结果。
    
> Step2: 计算远震的功率积分比值

1. 准备csv文件格式的远震目录，需要包含的字段有：
    - time，例如2006-06-20T10:02:07.78Z
    - latitude
    - longitude
    - depth(km)
    - mag 震级
    - dist(km) 震中距
    - f_min 目标频率范围的下限
    - f_max 目标频率范围的上限
    - Tb_B_time T_b时间窗开始的时间，例如2006-06-20T05:05:21.207727Z
    - Tb_E_time T_b时间窗结束的时间，例如2006-06-20T10:05:21.207727Z
    - Te_B_time T_e时间窗开始的时间，例如2006-06-20T10:07:10.520000Z
    - Te_E_time T_e时间窗结束的时间，例如2006-06-20T10:14:44.630000Z
   
   除了T_b和T_e相关的参数外，其他参数都比较容易获得。Test/bin/time_window.py可以用于一种确定T_b和T_e时间窗的方法，即T_b是远震到达前的某个固定时间长度，而T_e则是根据两个不同波速的震相的到时决定的。(该程序可同时输出震中距)
* 输入
    - T_b时间窗的长度 Tb
    - T_e开始时刻对应的震相的速度 vel_B
    - T_e结束时刻对应的震相的速度 vel_E
    - 台站的维度 sta_lat
    - 台站的经度 sta_lon
    - 已有的包含其他信息的远震目录 catalog_file
* 输出
    - 输入的远震目录文件增加了和时间窗相关的字段


2. 运行Test/bin/cal_tele_HighPIRation.py脚本，计算不同台站记录到的不同远震的功率积分比值
* 输入
    - 远震目录 teseismic_catalog
    - 功率积分数据库的路径 PIdatabase_path
    - 存放结果的路径 out_path
    - 台站列表 sta_list
    - 要计算的三分量列表 chn_list 

* 输出
    - 每个台站的结果存储为一个csv文件
    - 每个文件中每行为一个远震的计算结果，存储的字段有
        * event 远震发震时间
        * PI_b 远震到达前目标频率范围内的功率积分
        * PI_e 远震到达后目标频率范围内的功率积分
        * PIRatio 功率积分比PI_e/PI_b
        * PIRatio_log 功率积分比的对数值log10(PI_e/PI_b)

>> 下面的步骤是关于背景功率比值变化的计算和统计。由于计算背景功率积分比值时，应用的T_b和T_e时间窗是和远震相同的，所以只需将远震的发震时刻中的日期换成背景天的日期，然后生成相同的目录即可。
随后利用生成的目录计算功率谱积分比值即可。'

> Step3: 构造背景天目录

运行Test/bin/process_background_days/gen_background_catalog.py脚本，生成计算背景的功率积分比值时需要的目录
* 输入
    - 远震目录 teseismic_catalog
    - 保存输出结果的文件名称 out_file （xxx.csv）
    - 需要计算的背景的总天数 dayWindow
       
* 输出
    - 计算背景的功率积分比值的目录，与远震目录形式一致

> Step4: 计算背景的功率积分比值

运行Test/bin/process_background_days/cal_background_HighPIRation.py脚本，生成每个台站所有背景的功率积分比值。
*输入
    - step3中生成的背景目录 background_catalog
    - 功率积分数据库的路径 PIdatabase_path
    - 存放结果的路径 out_path
    - 台站列表 sta_list
    - 要计算的三分量列表 chn_list 

* 输出
    - 每个台站的结果存储为一个csv文件
    - 每个文件中每行为一个背景功率积分比值，存储的字段与远震的功率积分比值结果一致

> Step5: 归档远震与相应背景的功率积分比值

运行Test/bin/process_background_days/associate_background_PIR.py脚本，对于每个台站，将每个远震及其对应的背景功率积分比值数据集对应起来。
* 输入
    - 远震目录 teseismic_catalog
    - 背景功率积分比值路径 background_PIR_path
    - 输出路径 out_path
    - 需要计算的背景的总天数 dayWindow
    - 台站列表 sta_list
    - 要计算的三分量列表 chn_list 
    
* 输出
    - 每个台站的结果存储为一个csv文件
    - 文件中每行代表一个远震的结果，第一列为远震的发震时刻，后面所有列为该远震对应的背景功率积分比值结果。

> Step6: 计算置信水平 

运行Test/bin/cal_confidenceLevel.py脚本，每个远震在每个台站的数据中表现出来的动态触发的置信水平
* 输入
    - 存放step5结果的路径 background_PIR_associated_path
    - 存放远震功率积分比值的路径 tele_PIR_path
    - 输出结果的路径 out_path
    - 台站列表 sta_list
    - 要计算的三分量列表 chn_list 

* 输出
    - 每个台站的结果存储为一个csv文件
    - 文件中每行为一个远震的最终结果，存储的字段分别为：
        * tele 远震的发震时刻
        * confidence_level_value 置信水平值
        

#### utils
其他工具