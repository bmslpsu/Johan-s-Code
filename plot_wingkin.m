function plot_wingkin(raw_data)

    % Plot wing kinematics from tethered flight trackers:
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.theta_L),'r')
    hold off
    title('Left wing kinematics');
    ylabel('\theta')
    ylim([-30 30])
    subplot(3,1,2); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.eta_L1),'r')
    %plot(raw_data.frame_nr,radtodeg(raw_data.eta_L2),'b')
    hold off
    ylabel('\eta')
    ylim([0 90])
    subplot(3,1,3); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.phi_L),'r')
    hold off
    ylabel('\phi')
    ylim([-90 100])
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.theta_R),'b')
    hold off
    title('Right wing kinematics');
    ylabel('\theta')
    ylim([-30 30])
    subplot(3,1,2); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.eta_R1),'b')
    %plot(raw_data.frame_nr,radtodeg(raw_data.eta_R2),'r')
    hold off
    ylabel('\eta')
    ylim([0 90])
    subplot(3,1,3); hold on
    plot(raw_data.frame_nr,radtodeg(raw_data.phi_R),'b')
    hold off
    ylabel('\phi')
    ylim([-90 100])
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.theta_L(7500:8500)),'r')
    hold off
    title('Left wing kinematics, frame 7500 - 8500');
    ylabel('\theta')
    ylim([-30 30])
    subplot(3,1,2); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.eta_L1(7500:8500)),'r')
    %plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.eta_L2(7500:8500)),'b')
    hold off
    ylabel('\eta')
    ylim([0 90])
    subplot(3,1,3); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.phi_L(7500:8500)),'r')
    hold off
    ylabel('\phi')
    ylim([-90 100])
    hold off
    
    figure()
    hold on
    subplot(3,1,1); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.theta_R(7500:8500)),'b')
    hold off
    title('Right wing kinematics, frame 7500 - 8500');
    ylabel('\theta')
    ylim([-30 30])
    subplot(3,1,2); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.eta_R1(7500:8500)),'b')
    %plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.eta_R2(7500:8500)),'r')
    hold off
    ylabel('\eta')
    ylim([0 90])
    subplot(3,1,3); hold on
    plot(raw_data.frame_nr(7500:8500),radtodeg(raw_data.phi_R(7500:8500)),'b')
    hold off
    ylabel('\phi')
    ylim([-90 100])
    hold off

end

