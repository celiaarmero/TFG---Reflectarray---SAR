clear
clc
close all

%% Lectura de los datos de FEKO

DX = readmatrix('datos_feko_X.csv', ...
    'Delimiter',';','NumHeaderLines',1);

DY = readmatrix('datos_feko_Y.csv', ...
    'Delimiter',';','NumHeaderLines',1);

% Seleccion de la frecuencia de trabajo
DX = DX(abs(DX(:,3)-18e9)<1e6,:);
DY = DY(abs(DY(:,3)-18e9)<1e6,:);

%% Calculo de los coeficientes de reflexion complejos

RXX = DX(:,5).*exp(1i*deg2rad(DX(:,4)));
RYX = DX(:,8).*exp(1i*deg2rad(DX(:,7)));

RYY = DY(:,5).*exp(1i*deg2rad(DY(:,4)));
RXY = DY(:,8).*exp(1i*deg2rad(DY(:,7)));

%% Interpolacion de los coeficientes

FXX = interpolar(DX(:,1),DX(:,2),RXX);
FYX = interpolar(DX(:,1),DX(:,2),RYX);

FYY = interpolar(DY(:,1),DY(:,2),RYY);
FXY = interpolar(DY(:,1),DY(:,2),RXY);

%% Calculo de los coeficientes de los cinco paneles

for p = 1:5

    % Lectura de las dimensiones de los parches
    T = readtable(sprintf('parches_panel_%d_corregido.csv',p));

    % Lectura del archivo original de HADRA
    archivo = sprintf('Rho_%d_f18c00GHz.dat',p);

    rho = readmatrix(archivo,'CommentStyle','%');

    % Ordenar los elementos segun sus indices
    [~,orden] = ismember(rho(:,1:2),[T.NX,T.NY],'rows');

    % Calculo de los coeficientes de cada elemento
    for k = 1:size(rho,1)

        j = orden(k);

        if T.IFLAG(j) == 1

            a = T.a_final(j);
            b = T.b_final(j);

            XX = FXX(a,b);
            XY = FXY(a,b);
            YX = FYX(a,b);
            YY = FYY(a,b);

            rho(k,3:10) = [real(XX),imag(XX), ...
                          real(XY),imag(XY), ...
                          real(YX),imag(YX), ...
                          real(YY),imag(YY)];

        end

    end

    %% Lectura de la cabecera original

    fid = fopen(archivo,'r');
    cabecera = fgetl(fid);
    fclose(fid);

    %% Generacion del archivo de salida

    salida = sprintf('Rho_real_%d_f18c00GHz.dat',p);

    fid = fopen(salida,'w');

    fprintf(fid,'%s\n',cabecera);

    formato = ['%4d %4d',repmat(' %25.19f',1,8),'\n'];

    for k = 1:size(rho,1)

        fprintf(fid,formato,rho(k,:));

    end

    fclose(fid);

end

%% Funcion de interpolacion compleja

function F = interpolar(a,b,rho)

    Fre = scatteredInterpolant(a,b,real(rho),'linear','none');
    Fim = scatteredInterpolant(a,b,imag(rho),'linear','none');

    F = @(a,b) Fre(a,b) + 1i*Fim(a,b);

end