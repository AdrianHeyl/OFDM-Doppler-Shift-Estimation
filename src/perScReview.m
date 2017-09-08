%% perSC review
flag_perSC = 1;

if flag_perSC

    for k = 1:length(sc_active)
        sc_index = sc_active(k);

        figure;
        % accumulated phase rotation across OFDM symbols
        subplot(221);
        plot(deg_total(:,k)/pi);
        xlabel('symbol index');
        ylabel('phase rotation/pi');
        title(['phase info, N\_sc = ',num2str(N_sc),', Tsc #',num2str(sc_index)]);
        % unwraped phase
        subplot(222);
        plot(deg_unwraped_total(:,k)/pi);
        xlabel('symbol index');
        ylabel('phase rotation/pi');
        title('unwraped phase');

        % angular speed per symbol
        subplot(223);
        hold on;
        plot(w_symbol_diff(:,k));
        plot(w_symbol_gra(:,k));
        legend('diff','gradient');
        xlabel('symbol index');
        ylabel('\omega /pi');
        title('angular speed');

        % plot estimated speed per symbol
        subplot(223);
        hold on;
        plot(w_symbol_diff(:,k));
        plot(w_symbol_gra(:,k));
        legend('diff','gradient');
        xlabel('symbol index');
        ylabel('\omega /pi');
        title('angular speed');

        % plot estimated speed per symbol
        subplot(224);
        hold on;
        plot(v_symbol_gra(:,k));
        plot(v_symbol_gra_step2(:,k));
        legend('gra-ori','gra-MA');
        title(['sc #',num2str(sc_active(k))]);
    end
end
