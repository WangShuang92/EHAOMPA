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

% The EHAOMPA Optimization Algorithm
function [Best_FF,Best_P,conv]=EHAOMPA(N,T,LB,UB,Dim,F_obj)

Best_P=zeros(1,Dim);    
Best_FF=inf;
conv = zeros(1,T);
X=initialization(N,Dim,UB,LB);  


Xnew=X; 
Ffun=zeros(N,1);    
Ffun_new=zeros(N,1); 

t=0;

stepsize=zeros(N,Dim);

Xmin=repmat(ones(1,Dim).*LB,N,1);    
Xmax=repmat(ones(1,Dim).*UB,N,1);    

P=0.5;
FADs=0.2;

%------------------- Archive Initialization --------------------

archive_size=15;
archive_best_size=5;
archive_Positions = zeros(archive_size,Dim);
archive_fitness = zeros(1,archive_size);
archive_best_5position = zeros(archive_best_size,Dim);
archive_best_5fitness = zeros(1,archive_best_size);
sorted_positions = zeros(1,Dim);

%---------------------------------------------------------------

Exponent = 2;
sigma_initial=2;
sigma_final=0;
X_RH = zeros(1,Dim);


while t<T
    for i=1:size(X,1)   
        
        F_UB=X(i,:)>UB;
        F_LB=X(i,:)<LB;
        X(i,:)=(X(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
        
        
        Ffun(i)=F_obj(X(i,:));
        if Ffun(i)<Best_FF
            Best_FF=Ffun(i);
            Best_P=X(i,:);
        end
            
    end
    
    
    %------------------- Marine Memory saving -------------------
    
    if t==0
        Ffun_old=Ffun;    X_old=X;
    end
    
    Inx=(Ffun_old<Ffun);
    Indx=repmat(Inx,1,Dim);
    X=Indx.*X_old+~Indx.*X;
    Ffun=Inx.*Ffun_old+~Inx.*Ffun;
    
    Ffun_old=Ffun;    X_old=X;
    
    %------------------------------------------------------------    
    
    Elite = repmat(Best_P,N,1);

    CF=(1-t/T)^(2*t/T);
    
    to = 1:Dim;
    u = .0265;
    r0 = 10;
    r = r0 +u*to;
    omega = .005;
    phi0 = 3*pi/2;
    phi = -omega*to+phi0;
    x = r .* sin(phi);  % Eq. (6)
    y = r .* cos(phi); % Eq. (6)
    
    
    RL=0.05*Levy(N,Dim,1.5);
    
    
    %-------Select 15 positions from populations as archive ----------
    for j =1:archive_size
        archive_Positions(j,:) = X(j,:);   
    end
    
    %------------------- Calculate fitness -------------------
    for i =1:archive_size
        archive_fitness(1,i) = F_obj(archive_Positions(i,:));
    end
    
    %-------- Sort,select 5 best position from archive as Positions_R_Best------
    [sorted_archive_fitness,sorted_archive_indexes]=sort(archive_fitness);
    for newindex =1:archive_size
        sorted_positions(newindex,:) = archive_Positions(sorted_archive_indexes(newindex),:);
    end
    
    for i = 1:5
        archive_best_5position(i,:) = sorted_positions(i,:);
        archive_best_5fitness(1,i) = sorted_archive_fitness(i);
    end
    
    %------------------- Random populations -------------------
    rand_positions_index1 = floor(N*rand+1);
    rand_positions1 = X(rand_positions_index1,:);
    rand_positions_index2 = floor(N*rand+1);
    rand_positions2 = X(rand_positions_index2,:);
    
    %----------- Select a random pop from archive_best --------------
    rand_archiveBest_index = floor(archive_best_size*rand+1);
    archiveBest_rand = archive_best_5position(rand_archiveBest_index,:);
	
    %------------ Select a random pop from archive -------------------
    rand_archive_index = floor(archive_size*rand+1);
    archive_rand = archive_Positions(rand_archive_index,:);
    
    %------------------- Cauchy distribution -------------------
    for i = 1:N
        for j = 1:Dim
            cd = X(i,j) + 0.1*tan(pi*(rand-0.5));
        end
    end

    
    sigma = (((T-t)/(T-1))^Exponent) * (sigma_initial-sigma_final)+sigma_final;
    
  
    %-------------------------------------------------------------------------------------
    for i=1:size(X,1)
        %-------------------------------------------------------------------------------------

        if t<=(1/2)*T
		
            %------------------- Exploration -------------------

            if rand <0.5
                Xnew(i,:)=Best_P(1,:)*(1-t/T)+(mean(X(i,:))-Best_P(1,:))*rand(); % Eq. (22)
                X_RH =archiveBest_rand + cd*(X(i,:)-archive_rand)+sigma*(rand_positions1-rand_positions2); % Eq. (25)
                if F_obj(X_RH)<F_obj(Xnew(i,:))
                    Xnew(i,:) = X_RH;
                end
                Ffun_new(1,i)=F_obj(Xnew(i,:));
                if Ffun_new(i)<Ffun(i)
                    X(i,:)=Xnew(i,:);
                    Ffun(i)=Ffun_new(i);
                end
                
            else
                %-------------------------------------------------------------------------------------
                Xnew(i,:)=Best_P(1,:).*RL(i,:)+X((floor(N*rand()+1)),:)+(y-x)*rand;       % Eq. (24)
                X_RH =archiveBest_rand + cd*(X(i,:)-archive_rand)+sigma*(rand_positions1-rand_positions2); % Eq. (25)
                if F_obj(X_RH)<F_obj(Xnew(i,:))
                    Xnew(i,:) = X_RH;
                end
                Ffun_new(1,i)=F_obj(Xnew(i,:));
                if Ffun_new(i)<Ffun(i)
                    X(i,:)=Xnew(i,:);
                    Ffun(i)=Ffun_new(i);
                end
                
            end
            %-------------------------------------------------------------------------------------
            %--------------------- Exploitation ---------------------------------------
        else
            
            stepsize(i,:)=RL(i,:).*(RL(i,:).*Elite(i,:)-X(i,:));  % Eq. (27)
            Xnew(i,:)=Elite(i,:)+P*CF*stepsize(i,:);   % Eq. (26)
            Ffun_new(1,i)=F_obj(Xnew(i,:));
            if Ffun_new(i)<Ffun(i)
                X(i,:)=Xnew(i,:);
                Ffun(i)=Ffun_new(i);
            end
            Xnew_ROBL(i,:) = LB+UB-rand*X(i,:);   % Eq. (28)
            if F_obj(Xnew_ROBL(i,:)) < F_obj(X(i,:))
                X(i,:) = Xnew_ROBL(i,:);
            end
        end
     

    end
    
    
    %------------------ Detecting the best solution ------------------
    for i=1:size(X,1)   %
        
        F_UB=X(i,:)>UB;
        F_LB=X(i,:)<LB;
        X(i,:)=(X(i,:).*(~(F_UB+F_LB)))+UB.*F_UB+LB.*F_LB;
        
        
        Ffun(i)=F_obj(X(i,:));
        if Ffun(i)<Best_FF
            Best_FF=Ffun(i);
            Best_P=X(i,:);
        end
            
    end
    
    %---------------------- Marine Memory saving ----------------
    
    if t==0
        Ffun_old=Ffun;    X_old=X;
    end
    
    Inx=(Ffun_old<Ffun);
    Indx=repmat(Inx,1,Dim);
    X=Indx.*X_old+~Indx.*X;
    Ffun=Inx.*Ffun_old+~Inx.*Ffun;
    
    Ffun_old=Ffun;    X_old=X;
    

    %---------- Eddy formation and FADs' effect (Eq. 29) -----------
    
    if rand()<FADs
        U=rand(N,Dim)<FADs;
        X=X+CF*((Xmin+rand(N,Dim).*(Xmax-Xmin)).*U);
        
    else
        r=rand();  Rs=size(X,1);
        stepsize=(FADs*(1-r)+r)*(X(randperm(Rs),:)-X(randperm(Rs),:));
        X=X+stepsize;
    end
    
    
    t=t+1;
    conv(t)=Best_FF;
    
end

end


function [z] = Levy(n,m,beta)

    num = gamma(1+beta)*sin(pi*beta/2); 
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); 
    sigma_u = (num/den)^(1/beta);
    u = random('Normal',0,sigma_u,n,m);     
    v = random('Normal',0,1,n,m);
    z =u./(abs(v).^(1/beta));
  
  end


