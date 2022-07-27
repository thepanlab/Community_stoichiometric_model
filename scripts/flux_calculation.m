clear; clc
% copy and paste data from "mMole" tab in "balance_calculation"

data_table = readtable("dataset.csv")
data = data_table{:,:}
i=0;
j=0;
col=colormap(hsv(13))
t={'CD','CD+AM','CD+HM','CD+SR','CD+SR+SO4','CD+AM+HM','CD+AM+HM+SO4','CD+AM+SR','CD+AM+SR+SO4','CD+HM+SR','CD+HM+SR+SO4','CD+AM+HM+SR','CD+AM+HM+SR+SO4'}
% define overall reaction stoichiometries (lowest whole intergers on a mole basis)
a_table = readtable ("matrix.csv")
a = a_table{:,:}
 for c=1:13
  j=j+1;
  % set the present metabolisms that are being fit based on the populations present
  if c==1;     b=[a(:,1:3) zeros(8,4) a(:,8)];
  elseif c==2; b=[a(:,1:3) zeros(8,3) a(:,7:8)];
  elseif c==3; b=[a(:,1:3) zeros(8,2) a(:,6) zeros(8,1) a(:,8)];
  elseif c==4; b=[a(:,1:4) zeros(8,3) a(:,8)];
  elseif c==5; b=[a(:,1:5) zeros(8,2) a(:,8)];
  elseif c==6; b=[a(:,1:3) zeros(8,2) a(:,6:8)];
  elseif c==7; b=[a(:,1:3) zeros(8,2) a(:,6:8)];
  elseif c==8; b=[a(:,1:4) zeros(8,2) a(:,7:8)];
  elseif c==9; b=[a(:,1:5) zeros(8,1) a(:,7:8)];
  elseif c==10;b=[a(:,1:4) zeros(8,1) a(:,6) zeros(8,1) a(:,8)];
  elseif c==11;b=[a(:,1:6) zeros(8,1) a(:,8)];
  elseif c==12;b=[a(:,1:4) zeros(8,1) a(:,6:8)]; 
  elseif c==13;b=a(:,1:8); end
  for rep=1:3
    for d=1:7
      i=i+1;
      p(i,:)=data(((c-1)*7+d),[1 16 19 22 25 28 31 34]+rep-1); 
      p(i,5)=data(((c-1)*7+d),25+rep-1); 
      p(i,1)=-(data(((c-1)*7+d),1+rep-1)-data(((c-1)*7+d),4+rep-1)); 
      sol2(:,i) = lsqnonneg(b,p(i,:)');  cm(i,:)=b*sol2(1:end,i);
    end
    % PLOTS A plot the participating metabolisms
     subplot(3,5,j);plot([1:7],sol2(1,i-6:i),':r',[1:7],sol2(2,i-6:i),':b',[1:7],sol2(3,i-6:i),':g',[1:7],sol2(4,i-6:i),'--r',[1:7],sol2(5,i-6:i),'--b',[1:7],sol2(6,i-6:i),'--g',[1:7],sol2(7,i-6:i),'r',[1:7],sol2(8,i-6:i),'b','LineWidth',2); hold on;
    axis([0 7 0 1]);
    title(t{j})
  end
if c==1;
  % PLOTS A
  legend({'Homolactic acid','Eth/Ace','H2/CO2','Hydrogenic DVH','Sulfidogenic DVH','H2 based Sulfidogenic DVH','H2 based methanogen','Acetotrophic methanogen'});

end
end
