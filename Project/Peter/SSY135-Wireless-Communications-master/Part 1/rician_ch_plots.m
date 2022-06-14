function figure = rician_ch_plots( C,M,N )
%=====================================
colormap('jet')
% subplot(1, 2, 1)
% mesh(1:M, 0:N-1, abs(C));
% ylabel('frequency [1/(NTs)]')
% xlabel('time [Ts]')
% zlabel('|C[k, n]|')
set(gca, 'Ylim', [0 N-1], 'Xlim', [1 M])
% subplot(1, 2, 2)
figure = surf(1:M, 0:N-1, abs(C), 'MeshStyle', 'row');
ylabel('frequency [1/(NTs)]')
xlabel('time [Ts]')
set(gca, 'Ylim', [0 N-1], 'Xlim', [0 M-1])
view(2)
end