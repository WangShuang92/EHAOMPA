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

function X=initialization(N,Dim,UB,LB)

B_no= size(UB,2); % numnber of boundaries

if B_no==1
    X=rand(N,Dim).*(UB-LB)+LB;
end

% If each variable has a different lb and ub
if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
end