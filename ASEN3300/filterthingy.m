freqz(filt1.tf.num,filt1.tf.den,1000,120e3)

y = filter(filt1.tf.num,filt1.tf.den,x);