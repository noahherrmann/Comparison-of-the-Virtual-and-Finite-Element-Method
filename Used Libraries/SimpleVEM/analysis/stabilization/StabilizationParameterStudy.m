clc
clear
close all

Foulder = 'Stabilization_Study';

% Parametri di processo per infittimento mesh
nnBeam = [4,2];
nnHoledPlate = [3,3,3];
nnCookMembrane = [3,3];

% Tau di tentativo iniziale
tau0 = 1;

% Opzioni per la funzione (opzionale)
opzioni = optimoptions('fminunc', 'Display', 'iter');

% Ciclo per infittimento progressivo della mesh
niter = 5;
for i = 1:niter
    if i == 1 
        % Definzione delle funzioni da cui ricavare il valore ottimo di tau
        funzione1 = @(tau) CBConcComparisonUU(nnBeam,tau); 
        funzione2 = @(tau) CBDistrComparisonUU(nnBeam,tau); 
        funzione3 = @(tau) SBComparisonUU(nnBeam,tau); 
        funzione4 = @(tau) HPComparisonUU(nnHoledPlate,tau); 
        funzione5 = @(tau) CMComparisonUU(nnCookMembrane,tau); 
    
        % Ricerca del minimo per ognuno dei problemi
        [bestTau(1,i), minError(1,i)] = fminunc(funzione1, tau0, opzioni);
        [bestTau(2,i), minError(2,i)] = fminunc(funzione2, tau0, opzioni);
        [bestTau(3,i), minError(3,i)] = fminunc(funzione3, tau0, opzioni);
        [bestTau(4,i), minError(4,i)] = fminunc(funzione4, tau0, opzioni);
        [bestTau(5,i), minError(5,i)] = fminunc(funzione5, tau0, opzioni);
    else
        % Definzione delle funzioni da cui ricavare il valore ottimo di tau
        funzione1 = @(tau) CBConcComparisonUU(nnBeam,tau); 
        funzione2 = @(tau) CBDistrComparisonUU(nnBeam,tau); 
        funzione3 = @(tau) SBComparisonUU(nnBeam,tau); 
        funzione4 = @(tau) HPComparisonUU(nnHoledPlate,tau); 
        funzione5 = @(tau) CMComparisonUU(nnCookMembrane,tau); 
    
        % Ricerca del minimo per ognuno dei problemi
        [bestTau(1,i), minError(1,i)] = fminunc(funzione1, bestTau(1,i-1), opzioni);
        [bestTau(2,i), minError(2,i)] = fminunc(funzione2, bestTau(2,i-1), opzioni);
        [bestTau(3,i), minError(3,i)] = fminunc(funzione3, bestTau(3,i-1), opzioni);
        [bestTau(4,i), minError(4,i)] = fminunc(funzione4, bestTau(4,i-1), opzioni);
        [bestTau(5,i), minError(5,i)] = fminunc(funzione5, bestTau(5,i-1), opzioni);
    end
    
    nn(1,i) = prod(nnBeam);
    nn(2,i) = prod(nnHoledPlate);
    nn(3,i) = prod(nnCookMembrane);
    nnBeam = 2*nnBeam;
    nnHoledPlate = nnHoledPlate + 2;
    nnCookMembrane = 2*nnCookMembrane;
end

save('bestTau', 'bestTau');
save('Error',"minError")

figure()
plot(nn(1,:),bestTau(1,:),'LineWidth',1.5);
title("Mensola con carico concentrato","FontSize",20,"Interpreter",'latex');
% saveImage(Folder,'CBConc')

figure()
plot(nn(1,:),bestTau(2,:),'LineWidth',1.5);
title("Mensola con carico distribuito","FontSize",20,"Interpreter",'latex');
% saveImage(Foulder,'CBDistr')

figure()
plot(nn(1,:),bestTau(3,:),'LineWidth',1.5);
title("Trave appoggiata con carico distribuito","FontSize",20,"Interpreter",'latex');
% saveImage(Foulder,'SB')

figure()
plot(nn(2,:),bestTau(4,:),'LineWidth',1.5);
title("Piastra forata","FontSize",20,"Interpreter",'latex');
% saveImage(Foulder,'HP')

figure()
plot(nn(3,:),bestTau(5,:),'LineWidth',1.5);
title("Membrana di Cook","FontSize",20,"Interpreter",'latex');
% saveImage(Foulder,'CM')