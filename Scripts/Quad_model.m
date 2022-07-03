clear; clc
% copy and paste data from "mMole" tab in "balance_calculation"
data=[
0.055335222	0.055038444	0.053851333	0.047296444	0.047207778	0.047385111	0.695465505	0.690251251	0.690251251	0	0	0	0	0	0	0.049658895	0.04165266	0.035537563	0	0	0	0	0	0	0.077906828	0.078629519	0.07162475	0.095564152	0.088776281	0.095628552	0.168919661	0.165497347	0.167317458	0.00014	0.00014	0.00014
0.092623967	0.083690956	0.087430356	0.137323422	0.131817222	0.135115622	1.112049574	1.116742402	1.114424956	0.003970529	0.004107923	0.00451576	0	0	0	0.131031783	0.122582561	0.112615687	0.042325921	0.045155618	0.044199802	0	0	0	0.26540415	0.238118458	0.250869846	0.221357979	0.240682812	0.229220091	0.27928771	0.275699964	0.279556209	0.00014	0.00014	0.00014
0.240225011	0.343474	0.303082544	0.183123467	0.185313533	0.183974667	1.271918577	1.282347083	1.268210663	0.004443241	0.003026044	0.004277389	0.003607371	0.003527026	0.003567198	0.262850618	0.264488602	0.265074395	0.079675917	0.075017456	0.076837571	0	0	0	0.508330262	0.504506429	0.512208898	0.431437343	0.43660153	0.440491031	0.649823989	0.654901116	0.659896375	0.00014	0.00014	0.00014
0.733512811	0.783965033	0.743959389	0.208680911	0.208663178	0.209549844	1.382889472	1.385670407	1.381730749	0.003540658	0.005640123	0.005654907	0.00535432	0.005631117	0.005607981	0.31269811	0.378626584	0.325079049	0.138483931	0.14259898	0.138265051	0.003538154	0.003514395	0.003648261	0.698369253	0.653634784	0.69177853	0.583419197	0.617170723	0.606474546	0.765740924	0.759111362	0.781204707	0.00014	0.00014	0.00014
0.899039233	0.901057322	0.900760544	0.086626089	0.089738289	0.088656556	1.542642603	1.544438624	1.59530656	0.00302819	0.003099273	0.00330188	0.002910634	0.002985112	0.003086516	0.400529754	0.402025054	0.400735946	0.16280401	0.166132325	0.150683125	0.004077842	0.004002876	0.004315613	0.950966435	0.991279984	1.009941498	0.933451812	0.978379689	0.935495886	0.970275685	0.981147891	1.000845756	0.00014	0.00014	0.00014
1.015686356	1.026934233	1.027052944	0.088500667	0.088846467	0.087365733	1.703322712	1.703438584	1.725396384	0.002215284	0.002181264	0.002381034	0.002127713	0.002164971	0.003306982	0.578301879	0.562182655	0.567880059	0.141911404	0.135795183	0.131858141	0.004719978	0.004363175	0.004622721	1.108495884	1.125276534	1.134978112	1.167969184	1.163654508	1.124635606	1.160880887	1.121194449	1.199954299	0.00014	0.00014	0.00014
1.1544731	1.168273267	1.177681122	0.091545644	0.091829378	0.087351711	1.978079094	2.011102698	1.996676597	0.001587678	0.001593796	0.001787169	0.00146626	0.001809256	0.00180967	0.560767888	0.540567364	0.541293109	0.133748436	0.127128121	0.13570408	0.012936623	0.01277368	0.01265799	1.292164204	1.290458904	1.318065463	1.277414376	1.27934784	1.291145633	1.361856711	1.35630781	1.403454108	0.00014	0.00014	0.00014
																																			
0.038715667	0.033967222	0.044651222	0.052084444	0.052084444	0.051197778	0.70937018	0.713715391	0.731385916	0	0	0	0	0	0	0.042656665	0.044296824	0.042744569	0	0	0	0	0	0	0.073792541	0.076530753	0.076634316	0.086829266	0.088978224	0.089722264	0.166817866	0.165850179	0.169123928	0.00014	0.00014	0.00014
0.055348678	0.0673385	0.074639233	0.133368889	0.124413556	0.139841556	0.96280606	0.980042063	0.974856778	0.004445909	0.004325402	0.004482579	0	0	0	0.130432107	0.117835043	0.126480476	0.044426862	0.043810375	0.04357344	0	0	0	0.220275772	0.247469827	0.207138832	0.211881725	0.223084042	0.209138567	0.224462892	0.264953902	0.264002323	0.180753371	0.182555313	0.177014511
0.288421722	0.333116456	0.2665492	0.1936216	0.189090733	0.1887804	1.089072098	1.120154841	1.117750491	0.003386038	0.00335202	0.003396648	0.000472861	0.000572033	0.000522447	0.233314454	0.202485991	0.236182837	0.076215069	0.080753676	0.077755764	0.009970446	0.01284249	0.009838203	0.429463979	0.430358534	0.433719704	0.45612543	0.470718434	0.441048079	0.521573284	0.518354884	0.509684399	0.378120861	0.37982157	0.377850907
0.770728744	0.773221678	0.761261533	0.215171311	0.212227578	0.211739911	1.219625411	1.229851141	1.219567475	0.002947792	0.002928104	0.003101041	0.001950591	0.00200571	0.001803607	0.272914799	0.276140616	0.290086271	0.144953153	0.1519577	0.153824106	0.004123837	0.003272381	0.00098382	0.550922101	0.535475404	0.537364935	0.52780083	0.529500699	0.514431816	0.755157914	0.755166951	0.755463448	0.57642419	0.577389275	0.572885546
0.882686778	0.893519167	0.874050544	0.091352022	0.090367822	0.091653489	1.370746056	1.37338215	1.380537265	0.002054481	0.001967219	0.002075485	0.000929836	0.001023475	0.000981672	0.382219966	0.365349881	0.362743048	0.169064587	0.156247201	0.146552367	0.004517739	0.004454151	0.003874558	0.678989885	0.684702828	0.663163751	0.893398481	0.900953111	0.897713575	0.951056593	0.972143299	0.964676154	0.77504584	0.778798197	0.771294607
1.018980589	1.030703311	1.009869511	0.0891036	0.086213067	0.091063133	1.522098445	1.527341667	1.538146758	0.001815209	0.001774341	0.001790829	0.001843531	0.00187827	0.001845694	0.515646734	0.536078089	0.516913562	0.186977603	0.189536899	0.183258512	0.007192443	0.006764715	0.006711645	0.916472338	0.914580702	0.923175665	0.980676419	1.028065633	1.021411067	1.123591641	1.131093652	1.122970162	0.972823884	0.974305255	0.967946719
1.144560722	1.131087011	1.148359478	0.095916911	0.096954311	0.094551444	1.906064464	1.902037902	1.897663723	0.001143734	0.001119812	0.00118231	0.002363813	0.00243261	0.002208461	0.489472355	0.490604145	0.492109781	0.22665856	0.235758236	0.234835666	0.005878379	0.00518797	0.004268612	1.066553345	1.078349252	1.064510003	0.980613001	1.008371751	1.022359743	1.198567838	1.195022848	1.244491886	1.172267768	1.174716474	1.165709016];
t={'CD','CD+AM','CD+HM','CD+SR','CD+SR+SO4','CD+AM+HM','CD+AM+HM+SO4','CD+AM+SR','CD+AM+SR+SO4','CD+HM+SR','CD+HM+SR+SO4','CD+AM+HM+SR','CD+AM+HM+SR+SO4'};
% define overall reaction stoichiometries (lowest whole intergers on a mole basis)
%   HL,EA,HF,HG,SG,HS,HM,AM,GS  Legend: HL is homolactate, HF is Acetogenic, HG is Hydrogenic, SG is Sulfidogenic, HS is Hydrogen-based sulfidogenic, HM is Hydrogen-based methanogen, AM is Acetotrophic methanogen, and GS is glucose sink
%   1  2  3   4   5  6  7
a=[-1,-1, 0,  0,  0, 0, 0;  % glucose
    0, 2, 1,  1,  0, 0,-1;  % acetate
    2, 0,-1, -1,  0, 0, 0;  % lactate
    0, 4, 2,  0, -2,-4, 0;  % hydrogen
%     0, 2, 1,  1,  0,-1, 1;  % Carbon dioxide
    0, 0, 0,  0,  0, 1, 1;  % Methane
    0, 0, 0,0.5,0.5, 0, 0]; % Sulfide

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
%                 if k==1;     sol2(:,i) = lsqnonneg([zeros(size(b,1),1) b(:,2:end)],p(i,:)');  
%                 elseif k==2; sol2(:,i) = lsqnonneg([b(:,1) zeros(size(b,1),1) b(:,3:end)],p(i,:)');  
%                 elseif k==3; sol2(:,i) = lsqnonneg([b(:,1:2) zeros(size(b,1),1) b(:,4:end)],p(i,:)');  
%                 elseif k==4; sol2(:,i) = lsqnonneg([b(:,1:3) zeros(size(b,1),1) b(:,5:end)],p(i,:)');  
%                 elseif k==5; sol2(:,i) = lsqnonneg([b(:,1:4) zeros(size(b,1),1) b(:,6:end)],p(i,:)');  
%                 elseif k==6; sol2(:,i) = lsqnonneg([b(:,1:5) zeros(size(b,1),1) b(:,7:end)],p(i,:)');   
%                 elseif k==7; sol2(:,i) = lsqnonneg([b(:,1:6) zeros(size(b,1),1)],p(i,:)');  
%                 elseif k==8; sol2(:,i) = lsqnonneg([b(:,1:7)],p(i,:)');  end
%                 cm(i,:)=b*sol2(1:end,i);
                if k==1;     [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1) zeros(size(b,1),1) b(:,3) zeros(size(b,1),1) b(:,5:end)],p(i,:)');  
                elseif k==2; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:2) zeros(size(b,1),1) b(:,4) zeros(size(b,1),1) b(:,6:end)],p(i,:)');  
                elseif k==3; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:2) zeros(size(b,1),2) b(:,5:end)],p(i,:)');  
                elseif k==4; [sol2(:,i,rep,k,c)] = lsqnonneg([b(:,1:6) zeros(size(b,1),1)],p(i,:)'); end % no methanogen control
%                 elseif k==4; [sol2(:,i,rep,k,c)] = lsqnonneg(b,p(i,:)'); end % unconstrained control
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