function stop = plotTrajectoryCallback(x, optimValues, state)
    % Callback per fmincon: aggiorna la traiettoria, controlli e costo a ogni iterazione.

    % Variabili globali condivise con il main
    global nx nu N h xc yc r FigTraj J_history

    stop = false; % non fermare l'ottimizzazione

    switch state
        case 'init'
            % Inizializza figura e memoria del costo
            FigTraj = figure('Name', 'Evoluzione ottimizzazione', ...
                             'NumberTitle', 'off', 'Color', 'w');
            clf;
            tiledlayout(FigTraj, 3, 1, 'TileSpacing','compact', 'Padding','compact');
            J_history = []; % resetta lo storico del costo

        case 'iter'
            % === Aggiorna storico del costo ===
            if isfield(optimValues, 'fval')
                J_history(end+1) = optimValues.fval; %#ok<AGROW>
            end

            % --- Ricostruisci stati e controlli da x (=z)
            z = x(:);
            X = zeros(nx, N+1);
            U = zeros(nu, N);
            for k = 1:N+1
                X(:,k) = z((k-1)*(nx+nu) + (1:nx));
            end
            for k = 1:N
                U(:,k) = z((k-1)*(nx+nu) + nx + (1:nu));
            end
            t = 0:h:N*h;

            % --- Aggiorna la figura
            figure(FigTraj); clf;
            tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

            %% === 1) Traiettoria XY ===
            nexttile(1);
            hold on; grid on; axis equal;
            plot(X(1,:), X(2,:), 'b-', 'LineWidth', 1.5);
            plot(X(1,1), X(2,1), 'go', 'MarkerFaceColor','g');
            plot(X(1,end), X(2,end), 'ro', 'MarkerFaceColor','r');

            % Ostacolo
            theta = linspace(0,2*pi,100);
            fill(xc + r*cos(theta), yc + r*sin(theta), [1 0.8 0.8], ...
                 'EdgeColor','r','LineWidth',1.2);

            xlabel('$x$ [m]','Interpreter','latex');
            ylabel('$y$ [m]','Interpreter','latex');
            title(sprintf('Traiettoria (iter = %d)', optimValues.iteration), ...
                  'Interpreter','latex');
            legend({'Traiettoria','$x_{start}$','$x_{final}$','Ostacolo'}, ...
                   'Interpreter','latex','Location','best');
            set(gca,'TickLabelInterpreter','latex','FontSize',14);

            %% === 2) Controlli ===
            nexttile(2);
            plot(t(1:end-1), U(1,:), 'r-', 'LineWidth', 1.5);
            hold on;
            plot(t(1:end-1), U(2,:), 'b--', 'LineWidth', 1.5);
            grid on;
            xlabel('$t$ [s]', 'Interpreter','latex');
            ylabel('$u_1, u_2$', 'Interpreter','latex');
            title('Controlli', 'Interpreter','latex');
            legend({'$u_1$', '$u_2$'}, 'Interpreter','latex', 'Location','best');
            set(gca, 'TickLabelInterpreter','latex', 'FontSize', 14);

            %% === 3) Storia del costo J ===
            nexttile(3);
            semilogy(1:length(J_history), J_history, 'k-o', ...
                     'MarkerFaceColor','k','LineWidth',1.2);
            grid on;
            xlabel('Iterazione', 'Interpreter','latex');
            ylabel('$J(z)$', 'Interpreter','latex');
            title('Evoluzione del costo', 'Interpreter','latex');
            set(gca,'TickLabelInterpreter','latex','FontSize',14);

            drawnow;

        case 'done'
            % Plot finale più pulito (facoltativo)
            figure(FigTraj);
            sgtitle('Ottimizzazione completata', 'Interpreter','latex');
    end
end
