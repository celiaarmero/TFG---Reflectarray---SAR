clear
clc
close all

%% Datos de FEKO

DX = readtable('datos_feko_X1.csv','Delimiter',';');
DY = readtable('datos_feko_Y1.csv','Delimiter',';');

% Seleccion de la frecuencia de trabajo
DX = DX(abs(DX.Frequency-18e9)<=1e6,:);
DY = DY(abs(DY.Frequency-18e9)<=1e6,:);

DX = sortrows(DX,{'a','b'});
DY = sortrows(DY,{'a','b'});

% Dimensiones en mm
a = DX.a*1000;
b = DX.b*1000;

% Fases en forma compleja
zX = exp(1i*deg2rad(DX.Phase_deg));
zY = exp(1i*deg2rad(DY.Phase_deg));

%% Interpolacion de los datos

FXre = scatteredInterpolant(a,b,real(zX),'linear','none');
FXim = scatteredInterpolant(a,b,imag(zX),'linear','none');

FYre = scatteredInterpolant(a,b,real(zY),'linear','none');
FYim = scatteredInterpolant(a,b,imag(zY),'linear','none');

FX = @(a,b) FXre(a,b) + 1i*FXim(a,b);
FY = @(a,b) FYre(a,b) + 1i*FYim(a,b);

% Limites de las dimensiones
limInf = [min(a),min(b)];
limSup = [max(a),max(b)];

opts = optimset('Display','off', ...
    'TolX',1e-7,'TolFun',1e-9, ...
    'MaxFunEvals',240,'MaxIter',160);

% Guardar soluciones repetidas
cache = containers.Map('KeyType','char','ValueType','any');

%% Matrices para representar los cinco paneles

A = [];
B = [];
EX_total = [];
EY_total = [];
ET_total = [];

%% Calculo de las dimensiones de los parches

for p = 1:5

    datos = readmatrix(sprintf('faseob_%d.dat',p), ...
        'CommentStyle','%');

    ELEM = datos(:,1);
    NX = datos(:,2);
    NY = datos(:,3);
    X = datos(:,4);
    Y = datos(:,5);

    faseX_ideal = datos(:,6);
    faseY_ideal = datos(:,7);
    IFLAG = datos(:,8);

    N = length(ELEM);

    a_final = NaN(N,1);
    b_final = NaN(N,1);

    faseX_real = NaN(N,1);
    faseY_real = NaN(N,1);

    errorX = NaN(N,1);
    errorY = NaN(N,1);
    errorTotal = NaN(N,1);

    for k = 1:N

        if IFLAG(k) == 1

            phiX = mod(faseX_ideal(k)+180,360)-180;
            phiY = mod(faseY_ideal(k)+180,360)-180;

            clave = sprintf('%.8f_%.8f',phiX,phiY);

            if isKey(cache,clave)

                ab = cache(clave);

            else

                % Diferencias de fase
                EX = rad2deg(angle(zX.*exp(-1i*deg2rad(phiX))));
                EY = rad2deg(angle(zY.*exp(-1i*deg2rad(phiY))));

                % Buscar las mejores combinaciones iniciales
                [~,orden] = sort(EX.^2 + EY.^2);

                mejor = inf;
                ab = [a(orden(1)),b(orden(1))];

                % Optimizacion de las dimensiones
                for s = 1:min(3,length(orden))

                    inicio = [a(orden(s)),b(orden(s))];

                    objetivo = @(v) error_fase(v,phiX,phiY, ...
                        FX,FY,limInf,limSup);

                    v = fminsearch(objetivo,inicio,opts);

                    v = min(max(v,limInf),limSup);

                    coste = objetivo(v);

                    if coste < mejor
                        mejor = coste;
                        ab = v;
                    end

                    if mejor < 1e-8
                        break
                    end

                end

                cache(clave) = ab;

            end

            % Dimensiones finales en metros
            a_final(k) = ab(1)/1000;
            b_final(k) = ab(2)/1000;

            % Fases obtenidas
            zx = FX(ab(1),ab(2));
            zy = FY(ab(1),ab(2));

            faseX_real(k) = rad2deg(angle(zx));
            faseY_real(k) = rad2deg(angle(zy));

            % Errores de fase
            errorX(k) = rad2deg(angle(zx*exp(-1i*deg2rad(phiX))));
            errorY(k) = rad2deg(angle(zy*exp(-1i*deg2rad(phiY))));

            errorTotal(k) = hypot(errorX(k),errorY(k));

        end

    end

    %% Guardar los resultados

    a_mm = a_final*1000;
    b_mm = b_final*1000;

    resultados = table(ELEM,NX,NY,X,Y, ...
        faseX_ideal,faseY_ideal,a_final,b_final, ...
        a_mm,b_mm,faseX_real,faseY_real, ...
        errorX,errorY,errorTotal,IFLAG);

    writetable(resultados, ...
        sprintf('parches_panel_%d_corregido.csv',p));

    %% Distribucion de los elementos de cada panel

    mapaA = NaN(max(NY),max(NX));
    mapaB = NaN(max(NY),max(NX));

    mapaEX = NaN(max(NY),max(NX));
    mapaEY = NaN(max(NY),max(NX));
    mapaET = NaN(max(NY),max(NX));

    for k = 1:N

        mapaA(NY(k),NX(k)) = a_mm(k);
        mapaB(NY(k),NX(k)) = b_mm(k);

        mapaEX(NY(k),NX(k)) = errorX(k);
        mapaEY(NY(k),NX(k)) = errorY(k);
        mapaET(NY(k),NX(k)) = errorTotal(k);

    end

    % Union de los cinco paneles
    A = [A mapaA];
    B = [B mapaB];

    EX_total = [EX_total mapaEX];
    EY_total = [EY_total mapaEY];
    ET_total = [ET_total mapaET];

end

%% Representacion de la dimension a

figure
imagesc(A)
axis equal tight
axis xy
colormap(jet)
colorbar

xlabel('NX')
ylabel('NY')
title('Distribución de la dimensión a (mm)')

%% Representacion de la dimension b

figure
imagesc(B)
axis equal tight
axis xy
colormap(jet)
colorbar

xlabel('NX')
ylabel('NY')
title('Distribución de la dimensión b (mm)')

%% Error de fase en X

figure
imagesc(EX_total)
axis equal tight
axis xy
colormap(jet)
colorbar

xlabel('NX')
ylabel('NY')
title('Error de fase en X (°)')

%% Error de fase en Y

figure
imagesc(EY_total)
axis equal tight
axis xy
colormap(jet)
colorbar

xlabel('NX')
ylabel('NY')
title('Error de fase en Y (°)')

%% Error de fase combinado

figure
imagesc(ET_total)
axis equal tight
axis xy
colormap(jet)
colorbar

xlabel('NX')
ylabel('NY')
title('Error de fase combinado (°)')

%% Funcion para minimizar el error de fase

function coste = error_fase(v,phiX,phiY,FX,FY,limInf,limSup)

    fuera = max(limInf-v,0) + max(v-limSup,0);

    v = min(max(v,limInf),limSup);

    zx = FX(v(1),v(2));
    zy = FY(v(1),v(2));

    if ~isfinite(zx) || ~isfinite(zy) || ...
            abs(zx)<1e-10 || abs(zy)<1e-10

        coste = 1e9 + 1e6*sum(fuera.^2);
        return

    end

    ex = angle(zx*exp(-1i*deg2rad(phiX)));
    ey = angle(zy*exp(-1i*deg2rad(phiY)));

    coste = ex^2 + ey^2 + 1e4*sum(fuera.^2);

end