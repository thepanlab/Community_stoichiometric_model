clear; clc
% copy and paste data from "mMole" tab in "balance_calculation"
data_table = readtable("quad_dataset.csv");
data = data_table{:,:};
t={'CD','CD+AM','CD+HM','CD+SR','CD+SR+SO4','CD+AM+HM','CD+AM+HM+SO4','CD+AM+SR','CD+AM+SR+SO4','CD+HM+SR','CD+HM+SR+SO4','CD+AM+HM+SR','CD+AM+HM+SR+SO4'};
% define overall reaction stoichiometries (lowest whole intergers on a mole basis)
%   HL,EA,HF,HG,SG,HS,HM,AM,GS  Legend: HL is homolactate, HF is Acetogenic, HG is Hydrogenic, SG is Sulfidogenic, HS is Hydrogen-based sulfidogenic, HM is Hydrogen-based methanogen, AM is Acetotrophic methanogen, and GS is glucose sink
%   1  2  3   4   5  6  7

a_table = readtable ("matrix_1.csv")
a = a_table{:,:}

f1=figure(); f2=figure(); col=colormap(hsv(13)); 
for k=1:4; % run the entire list then remove a metabolism from the list and run again
    j=0;
    for c=1:2 % this is the number of cultures to look for
        if c==1;    b=[a(:,1:3) zeros(6,2) a(:,6:7)];
        elseif c==2;b=a(:,1:7); end
        j=j+1; %b=a;
        for rep=1:3 % this is the number of replicates to look for
            i=0;
            for d=1:7 % this is the number of days to look for
                i=i+1;
                p(i,:)=data(((c-1)*7+d),[1 16 19 25 31 34]+rep-1);
                p(i,1)=-(data(((c-1)*7+d),1+rep-1)-data(((c-1)*7+d),4+rep-1)); % calculate the glucose equivalent consumed
                % evaluate the collection of reactions removing one to test its importance
                if k==1;     [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1) zeros(size(b,1),1) b(:,3) zeros(size(b,1),1) b(:,5:end)],p(i,:)');  
                elseif k==2; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:2) zeros(size(b,1),1) b(:,4) zeros(size(b,1),1) b(:,6:end)],p(i,:)');  
                elseif k==3; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:2) zeros(size(b,1),2) b(:,5:end)],p(i,:)');  
                elseif k==4; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:6) zeros(size(b,1),1)],p(i,:)'); end % no methanogen control
%               elseif k==4; [sol2(:,i,rep,k,c)] = lsqnonneg(b,p(i,:)'); end % unconstrained control
                cm(i,:)=b*sol2(1:end,i,rep,k,c);
            end
            R(rep)=goodnessOfFit(p(7,:)',cm(7,:)','NRMSE')
            
            % PLOTS A plot the participating metabolisms
            figure(f1); subplot(4,2,j+(2*(k-1)));plot([1:7],sol2(1,i-6:i,rep,k,c),':r',[1:7],sol2(2,i-6:i,rep,k,c),':b',[1:7],sol2(3,i-6:i,rep,k,c),'--g',[1:7],sol2(4,i-6:i,rep,k,c),'--r',[1:7],sol2(5,i-6:i,rep,k,c),'--b',[1:7],sol2(6,i-6:i,rep,k,c),'r',...
                                                      [1:7],sol2(7,i-6:i,rep,k,c),'b','LineWidth',2); hold on;
            % PLOTS B plot the measured and fit metabolite profiles
            figure(f2); subplot(4,2,j+(2*(k-1)));plot(1:7,cm(i-6:i,1),':r',1:7,cm(i-6:i,2),':b',1:7,cm(i-6:i,3),':g',1:7,cm(i-6:i,4),'--k',1:7,cm(i-6:i,5),'--r',1:7,cm(i-6:i,6),'--b','LineWidth',2); 
                                        hold on; plot(1:7,p(i-6:i,1),'*r',1:7,p(i-6:i,2),'*b',1:7,p(i-6:i,3),'*g',1:7,p(i-6:i,4),'ok',1:7,p(i-6:i,5),'or',1:7,p(i-6:i,6),'ob')
        end
        % PLOTS A
        figure(f1); axis([0 7 0 2]); xlabel('Days'); ylabel('Flux (mmoles)'); title({t{c};['NRMSE ' num2str(round(mean(R),3)) ' +/- ' num2str(round(std(R),3))]})
        % PLOTS B
        figure(f2); axis([0 7 -1 1.5]); xlabel('Days'); ylabel({'Metabolite produced';'(mmoles)'}); title(t{c})       
        m(k,c)=mean(R)
        if c==12&k==4;
            % PLOTS A
            figure(f1); legend({'Homolactic acid','H2/CO2','Hydrogenic DVH','Sulfidogenic DVH','H2 based Sulfidogenic DVH','H2 based methanogen','Acetotrophic methanogen'});
            % PLOTS B
            figure(f2); legend({'Glucose','Acetate','Lactate','Hydrogen','Methane','Sulfide','Glucose','Acetate','Lactate','Hydrogen','Methane','Sulfide'});
        end
    end
    clear cm dml  
end
