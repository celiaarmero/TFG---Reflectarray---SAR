clear
clc
close all

%% Frecuencias de estudio

frecuencias = [17.75 18 18.25];
nombres = {'17c75','18c00','18c25'};

%% Bucle para cada frecuencia

for i = 1:length(frecuencias)

    % Lectura de los datos de FEKO
    magco = readmatrix(sprintf('magco%d.dat',i), ...
        'NumHeaderLines',2);

    faseco = readmatrix(sprintf('faseco%d.dat',i), ...
        'NumHeaderLines',2);

    a_FEKO = magco(:,1);

    % Coeficiente de reflexion complejo
    rho_co = 10.^(magco(:,2)/20).* ...
        exp(1i*deg2rad(faseco(:,2)));

    %% Bucle para los cinco paneles

    for p = 1:5

        % Lectura de las dimensiones de los parches
        datos = readmatrix( ...
            sprintf('tam_cuadrado_%d.dat',p));

        a = datos(:,6);

        % Lectura del archivo original de HADRA
        archivo = sprintf('Rho_%d_f18c00GHz.dat',p);

        rho = readmatrix(archivo,'CommentStyle','%');

        % Ordenar los elementos segun sus indices
        [~,orden] = ismember(rho(:,1:2), ...
            datos(:,2:3),'rows');

        a = a(orden);

        %% Calculo de los coeficientes de reflexion

        for k = 1:length(a)

            [~,idx] = min(abs(a_FEKO-a(k)));

            % Aproximacion folded reflectarray
            XX = rho_co(idx);
            YY = -XX;
            XY = 0;
            YX = 0;

            rho(k,3:10) = [real(XX),imag(XX), ...
                real(XY),imag(XY), ...
                real(YX),imag(YX), ...
                real(YY),imag(YY)];

        end

        %% Lectura de la cabecera original

        fid = fopen(archivo,'r');
        cabecera = fgetl(fid);
        fclose(fid);

        %% Exportacion de los resultados

        salida = sprintf( ...
            'Rho_folded_ideal_%d_f%sGHz.dat', ...
            p,nombres{i});

        fid = fopen(salida,'w');

        fprintf(fid,'%s\n',cabecera);

        formato = ['%4d %4d',repmat(' %25.19f',1,8),'\n'];

        for k = 1:size(rho,1)
            fprintf(fid,formato,rho(k,:));
        end

        fclose(fid);

    end

end