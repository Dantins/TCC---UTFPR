clear; clc;

%% Inicialização do hardware
run('tclab.m');  % Arquivo com as funções T1C(), h1(), led()
% Garante o desligamento seguro do TCLab mesmo em caso de erro ou interrupção (Ctrl+C)
cleanupObj = onCleanup(@() desligarTCLab());

%% Parâmetros de Simulação
run_time = 20;                % Tempo total da simulação (min)
loops = 60 * run_time;        % Número de ciclos (1 ciclo por segundo)
dt = 1.0;                     % Intervalo de tempo em segundos

% Inicialização das variáveis
tm = zeros(loops,1);          % Vetor para armazenar o tempo (s)
T1 = zeros(loops,1);          % Temperatura real (°C)
T1_ModeloSmith = zeros(loops,1);  % Temperatura do modelo Smith (°C)
Q1 = ones(loops,1) * 50;      % Potência fixa em 50%
% (u não é utilizado neste exemplo)

% Aplicar potência inicial ao hardware
h1(Q1(1));
Q = 50;  % Potência aplicada aos modelos

%% Criação do Gráfico
figure('Position',[100 100 800 600]);
drawnow;

%% Loop principal de aquisição e simulação
start_time = tic;
prev_time = toc(start_time);

try
    for i = 1:loops
        % Garantir que o tempo entre iterações seja aproximadamente dt segundos
        sleep_time = max(0.01, dt - (toc(start_time) - prev_time));
        pause(sleep_time);
        
        % Atualizar tempo
        t = toc(start_time);
        prev_time = t;
        tm(i) = t;
        
        % Ler temperatura real do TCLab
        T1(i) = T1C();
        
        % --- Modelo Smith ---
        if i == 1
            % Condição inicial
            T1_ModeloSmith(i) = T1(1);
        elseif i <= 25
            % Enquanto o atraso não ocorre (i-25 ≤ 0), mantém o valor anterior
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1);
        else
            % A atualização ocorre para i > 25.
            % Utiliza heaviside(i-25-eps) para garantir que o termo seja 0 quando i-25 é zero
            T1_ModeloSmith(i) = T1_ModeloSmith(i-1) + 0.00340 * exp(-tm(i-25) / 174) * Q * heaviside(i-25 - eps);
        end
        
        % --- Atualização Gráfica ---
        clf;
        subplot(2,1,1);
        plot(tm(1:i), T1(1:i), 'r-', 'LineWidth', 2); hold on;
        plot(tm(1:i), T1_ModeloSmith(1:i), 'g--', 'LineWidth', 2);
        ylabel('Temperatura (°C)', 'FontSize', 14);
        legend('Temperatura Real', 'Modelo Smith', 'FontSize', 12);
        grid on;
        
        subplot(2,1,2);
        plot(tm(1:i), Q1(1:i), 'r-', 'LineWidth', 2);
        ylabel('Potência (%)', 'FontSize', 14);
        xlabel('Tempo (s)', 'FontSize', 14);
        legend('Potência (%)', 'FontSize', 12);
        grid on;
        
        drawnow;
    end
    
    % --- Cálculo do Erro e Salvamento ---
    erro_medio_Smith = mean(abs(T1 - T1_ModeloSmith));
    fprintf('Erro Médio Absoluto entre Temperatura Real e Modelo Smith: %.2f °C\n', erro_medio_Smith);
    
    % Salvar gráfico em PNG
    saveas(gcf, 'Grafico_TCLab_Modelo_Smith.png');
    
    % Salvar dados em arquivo TXT
    dados = [tm, T1, T1_ModeloSmith, Q1];
    header = {'Tempo_s', 'Temperatura_Real_C', 'Modelo_Smith_C', 'Potencia_pct'};
    T_data = array2table(dados, 'VariableNames', header);
    writetable(T_data, 'Dados_simulacao_08_Smith.txt', 'Delimiter', '\t');
    
catch ME
    disp(['Erro: ' ME.message]);
    rethrow(ME);
end

disp('Simulação concluída.');

%% Função de desligamento seguro do TCLab
function desligarTCLab()
    try
        h1(0);    % Desliga o aquecedor
        led(0);   % Apaga o LED
        disp('Dispositivos TCLab desligados com segurança.');
    catch
        disp('Erro ao tentar desligar os dispositivos.');
    end
end
