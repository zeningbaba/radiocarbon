function [meanage,meanProjVentAge,meanProjDD14Ccorr,ProjVentAge_ellipse, ...
    ProjDD14Ccorr_ellipse] = Projection_Age(age_radiocarbon)

%calculate error ellipse for radiocarbon
% update 2022.8.1 written by Tianyu Chen (tianyuchen@nju.edu.cn)



marine20;

Range = 15000; % serch range
scalefactor = 5;
N_age=500;
N_14C=500;
cal=sortrows(cal,1);

% convert to projection 2 +400 year
% cal(:,2) = cal(:,2)+400;

% smooth data
% cal(:,2) = smoothdata(cal(:,2),'gaussian',60);
% 
% cal(:,4) = smoothdata(cal(:,4),'gaussian',60);
%
% calcualte decay corrected atmosphere radiocarbon age
 AtmR(:,1) = cal(:,1);
 AtmR(:,2) = cal(:,2)-5568/5730*cal(:,1);
 AtmR(:,3) = cal(:,3);
 
 AtmRDecay = NaN(Range/scalefactor,3);
AtmRDecay(:,1) =[scalefactor:scalefactor:Range]';
AtmRDecay(:,2) = interp1(AtmR(:,1),AtmR(:,2),AtmRDecay(:,1));
MCAtmRDecay = NaN(Range/scalefactor,N_age);
for tt = 1:Range/scalefactor
[c, index] = min(abs(AtmRDecay(tt,1)-AtmR(:,1)));
AtmRDecay(tt,3) = AtmR(index,3);
MCAtmRDecay(tt,:) = AtmRDecay(tt,2) + AtmRDecay(tt,3).*randn(N_age,1)';
end
%         end
    age_radiocarbon=sortrows(age_radiocarbon,1);
    smplnumber = size(age_radiocarbon,1);
    %     age_radiocarbon(:,2)=0;
    age_radiocarbon(:, 2) =  age_radiocarbon(:, 2)/2;
    age_radiocarbon(:, 4) =  age_radiocarbon(:, 4)/2;
    clear ProjVentAge_ellipse AgeRand RadioAgeRand
    
    AgeRand = NaN(N_age, smplnumber);
    RadioAgeRand = NaN(N_14C, smplnumber);
    
    
        for i = 1 : smplnumber
           
            AgeRand(:,i) = age_radiocarbon(i,1) + age_radiocarbon(i,2).*randn(N_age,1); % 1 sigma
            
            if AgeRand(:,i)<0
                warning('There are negative calendar ages in sample %d of %s during monte carlo simulation', i,loc)
                AgeRand(AgeRand(:,i)<0,i)=0;
                disp('Negative ages have been forced to equal zero in the ventilation calculation')
            end
            
            
            RadioAgeRand(:,i) = age_radiocarbon(i,3) + age_radiocarbon(i,4).*randn(N_14C,1); % 1 sigma
            
           
%            smplRDecayTrack = RadioAgeRand(:,i)-5568/5730*[scalefactor:scalefactor:Range]; %
           smplRDecayTrack = RadioAgeRand(:,i)-5568/5730*[scalefactor:scalefactor:Range]; 
            kk = NaN(N_14C,N_14C);
            ProjVentilation=NaN(N_14C,N_14C);
            ProjCalendar=NaN(N_14C,N_14C);
           for pp = 1:N_14C
           DifR = smplRDecayTrack'-MCAtmRDecay(:,pp);
%            EasyArray = DifR<0|DifR==0;
%            [Y,kk(pp,:)] = max(EasyArray,[],1);
           EasyArray = abs(DifR);
           [Y,kk(pp,:)] = min(EasyArray);
           clear EasyArray DifR  
%            ProjVentilation(pp,:)=kk(pp,:)*scalefactor-AgeRand(:,i)';
           ProjVentilation(pp,:)=kk(pp,:)*scalefactor-AgeRand(:,i)';
           ProjCalendar(pp,:)= AgeRand(:,i)';
           end
           
           clear smplRDecayTrack kk
                                
                             
            [meanage(i), meanProjVentAge(i),ellipse_ProjVentAge_x(:,i),ellipse_ProjVentAge_y(:,i)] = ...
                error_ellipse(reshape(ProjCalendar(:,:),[],1),reshape(ProjVentilation(:,:),[],1));
         clear ProjCalendar ProjVentilation kk
        end
       
        
     ProjVentAge_ellipse = NaN(size(ellipse_ProjVentAge_x,1),2*i);
     ProjVentAge_ellipse(:,1:2:end) = ellipse_ProjVentAge_x; 
     ProjVentAge_ellipse(:,2:2:end) = ellipse_ProjVentAge_y;
     
    % ProjVentAge_ellipse(:,2:2:end) = ellipse_ProjVentAge_y+400; %+400 for project2
     
     ProjDD14Ccorr_ellipse = NaN(size(ellipse_ProjVentAge_x,1),2*i);
     ProjDD14Ccorr_ellipse(:,1:2:end) = ellipse_ProjVentAge_x;
     ProjDD14Ccorr_ellipse(:,2:2:end) = (exp(-(ellipse_ProjVentAge_y)*log(2)/5568)-1)*1000; 
%      ProjDD14Ccorr_ellipse(:,2:2:end) = (exp(-(ellipse_ProjVentAge_y+400)*log(2)/5568)-1)*1000; %+400 for project2
%      meanProjDD14Ccorr = (exp(-(meanProjVentAge+400)*log(2)/5568)-1)*1000; %+400 for project2
     meanProjDD14Ccorr = (exp(-(meanProjVentAge)*log(2)/5568)-1)*1000;

end




