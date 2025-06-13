clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                 % Tempo total da simulação (min)
loops = 60 * run_time;         % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                      % Intervalo de tempo (s)

% Inicialização das variáveis
tm = zeros(loops,1);           % Vetor de tempo (s)
T1 = zeros(loops,1);           % Temperatura real (°C)
T1_ModeloNL = zeros(loops,1);  % Temperatura modelada não linear (não utilizada aqui)
T1_ModeloL = zeros(loops,1);   % Temperatura modelada linear (não utilizada aqui)
Q1 = ones(loops,1) * 50;       % Potência fixa em 50%
u = zeros(loops,1);            % Entrada no modelo linear (se necessário)

% Aplicar potência inicial ao hardware
h1(Q1(1));

%% Criação do Gráfico
figure('Position',[100 100 800 600]);
drawnow;

%% Loop principal
start_time = tic;
prev_time = toc(start_time);

try
    for i = 1:loops
        % Tempo de espera: garante que o intervalo seja de dt segundos
        sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
        pause(sleep_time);
        
        % Atualizar tempo
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        % Ler temperatura real do TCLab
        T1(i) = T1C();
        
        % Atualizar gráfico
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'FontSize', 12);
        grid on;
        
        subplot(2,1,2);
        plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
        ylabel('Potência (%)', 'FontSize', 14);
        xlabel('Tempo (s)', 'FontSize', 14);
        legend('Potência (%)', 'FontSize', 12);
        grid on;
        
        drawnow;
    end

    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Modelo_Malha_Aberta_50.png');

    % Salvar dados em arquivo TXT
    dados = [tm, T1, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Potencia_pct'};
    T = array2table(dados, 'VariableNames', header);
    writetable(T, 'Dados_simulacao_Malha_Aberta_50.txt', 'Delimiter', '\t');
    
catch ME
    disp(['Erro: ', ME.message]);
    rethrow(ME);
end

disp('Simulação concluída.');

%% Função de desligamento seguro
function desligarTCLab()
    try
        h1(0);    % Desliga o aquecedor
        led(0);   % Apaga o LED
        disp('Dispositivos TCLab desligados com segurança.');
    catch
        disp('Erro ao tentar desligar os dispositivos.');
    end
end
