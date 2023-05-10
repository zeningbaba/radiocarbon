function [meanage,smplnumber,meandelta14C,meandeltadelta14C,deltaellipse,deltadeltaellipse] = Delta14C( A,tag)
%UNTITLED2 Summary of this function goes here
% updated by Tao Li (taoli@nigpas.ac.cn) 2023.5.10
% written by Tianyu Chen (tianyuchen@nju.edu.cn) 
% 
marine20;
A=sortrows(A,1);
smplnumber = length(A(:,1));
A(:, 2) =  A(:, 2)/2;
A(:, 4) =  A(:, 4)/2;

N_age=1000;
N_14C=1000;
AgeRand = NaN(N_age, smplnumber);
delta14CATM = NaN(N_age, smplnumber);
sigmadeltaCATM = NaN(1,smplnumber);
FmodernRand = NaN(N_14C, smplnumber);
RadioAgeRand = NaN(N_14C, smplnumber);
MCdelta14Catm = NaN(N_14C,N_age,smplnumber);
MCdelta14Csmpl = NaN(N_14C,N_age,smplnumber);
MCdeltadelta14C = NaN(N_14C,N_age,smplnumber);


if strcmp(tag,'Y')
for i = 1 : smplnumber; 
   AgeRand(:,i) = A(i,1) + A(i,2).*randn(N_age,1); % 1 sigma
  FmodernRand(:,i) = A(i,3) + A(i,4).*randn(N_14C,1);% 1 sigma
    
  delta14CATM(:,i) = interp1(cal(:,1),cal(:,4),AgeRand(:,i));
  %deltadelta14Crand(:,i) = -delta14CrandATM(:,i)+delta14Crand(:,i);
  [c, index] = min(abs(cal(:,1)-A(i,1)));
  sigmadeltaCATM(i) = cal(index,5);
  
  
  %sigmadeltadeltaC(i) =  sqrt(sigmadeltaCATM(i)^2+(coal_A(i,4)/2)^2);
  for j = 1:N_age
  MCdelta14Csmpl(:,j,i) = (FmodernRand(:,i).*exp(AgeRand(j,i)/8267)-1)*1000;
  MCdelta14Catm(:,j,i) = delta14CATM(j,i) + sigmadeltaCATM(i).*randn(N_14C,1);
  MCdeltadelta14C(:,j,i) = MCdelta14Csmpl(:,j,i)-MCdelta14Catm(:,j,i);
  end
  
  [meanage(i), meandelta14C(i),ellipse_delta_x(:,i),ellipse_delta_y(:,i)] = ...
      error_ellipse(reshape(repmat(AgeRand(:,i),1,N_14C)',[],1),reshape(MCdelta14Csmpl(:,:,i),[],1));
  [meanage2(i), meandeltadelta14C(i),ellipse_deltadelta_x(:,i),ellipse_deltadelta_y(:,i)] = ...
      error_ellipse(reshape(repmat(AgeRand(:,i),1,N_14C)',[],1),reshape(MCdeltadelta14C(:,:,i),[],1));
 
end
 
else
    
for i = 1 : smplnumber; 
   AgeRand(:,i) = A(i,1) + A(i,2).*randn(N_age,1); % 1 sigma
 RadioAgeRand(:,i) = A(i,3) + A(i,4).*randn(N_14C,1); % 1 sigma
    
  delta14CATM(:,i) = interp1(cal(:,1),cal(:,4),AgeRand(:,i));
  %deltadelta14Crand(:,i) = -delta14CrandATM(:,i)+delta14Crand(:,i);
  [c, index] = min(abs(cal(:,1)-A(i,1)));
  sigmadeltaCATM(i) = cal(index,5);
  
  
  %sigmadeltadeltaC(i) =  sqrt(sigmadeltaCATM(i)^2+(coal_A(i,4)/2)^2);
  for j = 1:N_14C
  MCdelta14Csmpl(:,j,i) = (exp(-RadioAgeRand(:,i)/8033).*exp(AgeRand(j,i)/8267)-1)*1000;
  MCdelta14Catm(:,j,i) = delta14CATM(j,i) + sigmadeltaCATM(i).*randn(N_14C,1);
  MCdeltadelta14C(:,j,i) = MCdelta14Csmpl(:,j,i)-MCdelta14Catm(:,j,i);
  end
  
  [meanage(i), meandelta14C(i),ellipse_delta_x(:,i),ellipse_delta_y(:,i)] = ...
      error_ellipse(reshape(repmat(AgeRand(:,i),1,N_14C)',[],1),reshape(MCdelta14Csmpl(:,:,i),[],1));
  [meanage(i), meandeltadelta14C(i),ellipse_deltadelta_x(:,i),ellipse_deltadelta_y(:,i)] = ...
      error_ellipse(reshape(repmat(AgeRand(:,i),1,N_14C)',[],1),reshape(MCdeltadelta14C(:,:,i),[],1));
 
end

end

deltaellipse = NaN(size(ellipse_delta_x,1),2*size(ellipse_delta_x,2));
deltaellipse(:,1:2:end) = ellipse_delta_x;
deltaellipse(:,2:2:end) = ellipse_delta_y;
deltadeltaellipse = NaN(size(ellipse_deltadelta_x,1),2*size(ellipse_deltadelta_x,2));
deltadeltaellipse(:,1:2:end) = ellipse_deltadelta_x;
deltadeltaellipse(:,2:2:end) = ellipse_deltadelta_y;

clear AgeRand delta14Crand sigmadeltaC MCdelta14C AgeInput delta14Cinput deltadelta14Crand sigmadeltadeltaC MCdeltadelta14C AgeInput deltadelta14Cinput
clear MCdelta14CATM MCdelta14Csmpl FmodernRand
clear MCdelta14Catm RadioAgeRand delta14CATM

end
