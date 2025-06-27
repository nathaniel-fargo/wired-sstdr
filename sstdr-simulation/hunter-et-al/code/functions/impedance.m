function complex_Z = impedance(Q,f,aaa)

mag_info=aaa.data(1:201,2);mag_info=[1956;mag_info];
mag_info=mag_info(1:101);
phase_info=aaa.data(202:402,2);phase_info=[0;phase_info];
phase_info=phase_info(1:101);
freqs=aaa.data(1:201,1);freqs=[1;freqs];
freqs=freqs(1:101);

pos_mag_spline = interp1(freqs,mag_info,f(1:Q/2),'spline')';
neg_mag_spline=flipud(pos_mag_spline);
mag_spline=[pos_mag_spline;neg_mag_spline];
pos_phase_spline = interp1(freqs,phase_info,f(1:Q/2),'spline')';
neg_phase_spline=flipud(pos_phase_spline)*-1;
phase_spline = [pos_phase_spline;neg_phase_spline];
    
real_parts=mag_spline.*cos(phase_spline*(pi/180));real_parts(1,1)=1;
imag_parts=1j*mag_spline.*sin(phase_spline*(pi/180));
    
complex_Z=(real_parts+imag_parts);
end