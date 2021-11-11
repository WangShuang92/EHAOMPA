%___________________________________________________________________________________________________________%
%  Enhanced Hybrid Aquila Optimizer and Marine Predators Algorithm (EHAOMPA) source codes demo 1.0          %
%                                                                                                           %
%  Developed in MATLAB R2016a                                                                               %
%                                                                                                           %
%  Author and programmer: Shuang Wang, Heming Jia, Laith Abualigah, Guanjun Lin, Hongwei Wei, Qingxin Liu   %
%                                                                                                           %
%         e-Mail: wangshuang14@mails.ucas.ac.cn, wang_shuang9279@163.com                                    %
%                 jiaheming@fjsmu.edu.cn, jiaheminglucky99@126.com                                          %
%                                                                                                           %
%  Homepage: https://www.researchgate.net/profile/Shuang-Wang-52                                            %
%  Homepage: https://www.researchgate.net/profile/Heming-Jia                                                %
%___________________________________________________________________________________________________________%

clear all 
clc
close all
Solution_no=30;  
F_name='F1';    
M_Iter=500;    


[LB,UB,Dim,F_obj]=Get_Functions_details(F_name); 

[Best_FF,Best,conv]=EHAOMPA(Solution_no,M_Iter,LB,UB,Dim,F_obj);

 
 figure('Position',[454   445   694   297]);
 subplot(1,2,1);
 func_plot(F_name);
 title('Parameter space')
 xlabel('x_1');
 ylabel('x_2');
 zlabel([F_name,'( x_1 , x_2 )'])
 
 
 subplot(1,2,2);
 semilogy(conv,'Color','r','linewidth',1.5)


 title(F_name)
 xlabel('Iteration');
 ylabel('Best Score');
 axis tight
 legend('EHAOMPA');

 display(['The best solution obtained by EHAOMPA is  ', num2str(Best)]);
 display(['The best optimal values by EHAOMPA is : ', num2str(Best_FF)]);



