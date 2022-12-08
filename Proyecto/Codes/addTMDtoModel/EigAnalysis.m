function [T,Gphi] = EigAnalysis(wi,k,g)
% P.Heresi THAMDOF Function
%%%%%%%%%%%%%%

Nst = length(wi);

% M and K matrices
M = diag(wi)/g;         % Mass Matrix
k_aux = k(2:end);   k_aux(end+1) = 0;
K = diag(k+k_aux) - diag(k(2:end),1) - diag(k(2:end),-1);

[phi,w2] = eig(K,M);
w = sqrt(diag(w2));     % Undamped frequencies
[w,index] = sort(w);
T = 2*pi./w;            % Undamped periods

% Sort vectors (modal shapes) and normalize them at roof: phi_roof = 1.0
sphi = phi;
for i = 1:length(wi)
    sphi(:,i) = phi(:,index(i))/ phi(end,index(i));
end
phi = sphi;             % Normalized modal shapes

r = ones(Nst,1);
Gamma = phi'*M*r./diag(phi'*M*phi);

Gphi = phi;
for i = 1:Nst
    Gphi(:,i) = Gamma(i)*phi(:,i);
end


function [] = plotUndeformed(handles)
global Building g lengthT forceT

Nst = length(Building.Story);

plot(handles.BuildingAxes,[0;Building.do],[0;Building.H],'-o','linewidth',2)
axes(handles.BuildingAxes);
grid on
xlabel(['do [' lengthT ']'])
ylabel(['Height [' lengthT ']'])
ylim([0 1.2*Building.H(end)])

doAux = [0;Building.do];
HAux = [0;Building.H];
for i = 1:Nst
    Xtext = (doAux(i) + doAux(i+1))/2;
    Ytext = (HAux(i) + HAux(i+1))/2;
    text(Xtext,Ytext,['   K_' num2str(i) ' = ' num2str(Building.K(i)) ' ' forceT '/' lengthT])
    Xtext = Building.do(i);
    Ytext = Building.H(i);
    text(Xtext,Ytext,['   W_' num2str(i) ' = ' num2str(Building.W(i)) ' ' forceT])
end
