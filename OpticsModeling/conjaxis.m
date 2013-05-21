function b=conjaxis(a,padlen)
% F=CONJAXIS(T,PADLEN) returns the discrete conjugate axis F from given axis T. By
% By conjugate, we mean that both axes represent independent variables related
% by discrete Fourier transfrom. PADLEN is the ratio between lengths of F
% and T.
% May 2010, Shalin Mehta (shalin.mehta@gmail.com)

    da=a(2)-a(1); %We assume uniform sampling. Compute the sampling rate.
    bcut=1/(2*da);
    b=linspace(-bcut,bcut,round(padlen*length(a)));
    bzeroidx=floor(0.5*length(b))+1;
    b=b-b(bzeroidx); %This ensures that zero value occurs at the right position.
end