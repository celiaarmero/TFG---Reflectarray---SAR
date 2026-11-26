clear
clc
close all

%% Datos de FEKO a 18 GHz

FEKO = readmatrix('faseco2.dat','NumHeaderLines',2);

a_FEKO = FEKO(:,1);
fase_FEKO = FEKO(:,2);

%% Matrices para las figuras

A = [];
Error = [];

%% Calculo de las dimensiones de los cinco paneles

for p = 1:5

    datos = readmatrix(sprintf('faseob_%d.dat',p), ...
        'NumHeaderLines',1);

    ELEM = datos(:,1);
    NX = datos(:,2);
    NY = datos(:,3);
    X = datos(:,4);
    Y = datos(:,5);

    fase_ideal = datos(:,6);

    N = length(ELEM);

    a = zeros(N,1);
    fase_real = zeros(N,1);
    error_fase = zeros(N,1);

    % Busqueda del tamaño de cada parche
    for k = 1:N

        diferencia = abs(mod(fase_FEKO - ...
            fase_ideal(k)+180,360)-180);

        [error_fase(k),idx] = min(diferencia);

        a(k) = a_FEKO(idx);
        fase_real(k) = fase_FEKO(idx);

    end

    % Parches cuadrados
    b = a;

    %% Guardar los resultados

    resultados = [ELEM NX NY X Y a b ...
        fase_ideal fase_real error_fase];

    writematrix(resultados, ...
        sprintf('tam_cuadrado_%d.dat',p), ...
        'Delimiter','tab');

    %% Calculo de los errores

    fprintf('\nPanel %d\n',p);
    fprintf('Error medio: %.3f grados\n',mean(error_fase));
    fprintf('Error maximo: %.3f grados\n',max(error_fase));
    fprintf('Error RMS: %.3f grados\n',sqrt(mean(error_fase.^2)));

    %% Distribucion de los elementos

    mapaA = zeros(max(NY),max(NX));
    mapaError = zeros(max(NY),max(NX));

    for k = 1:N

        mapaA(NY(k),NX(k)) = a(k)*1000;
        mapaError(NY(k),NX(k)) = error_fase(k);

    end

    % Union de los cinco paneles
    A = [A mapaA];
    Error = [Error mapaError];

end

%% Representacion de las dimensiones

figure
imagesc(A)
axis equal tight
axis xy
colormap(jet)

clim([min(a_FEKO) max(a_FEKO)]*1000)

c = colorbar;
c.Label.String = 'a = b (mm)';

xlabel('NX')
ylabel('NY')
title('Distribución de las dimensiones a = b')

%% Representacion del error de fase

figure
imagesc(Error)
axis equal tight
axis xy
colormap(jet)

c = colorbar;
c.Label.String = 'Error de fase (°)';

xlabel('NX')
ylabel('NY')
title('Distribución del error de fase')