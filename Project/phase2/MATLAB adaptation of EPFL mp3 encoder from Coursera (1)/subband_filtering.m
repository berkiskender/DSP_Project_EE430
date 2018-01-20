function s = subband_filtering(x,h)
    d=zeros(512,1);
    cq=zeros(64,1);
    sk=zeros(32,1);
    for i=1:512
        d(i)=h(i).*x(i);
    end
    for q=1:64
        for p=1:8
        cq(q)=cq(q)+((-1)^(p-1))*d(64*(p-1)+q);
        end
    end
    for k=1:32
        for q=1:64
            sk(k)=sk(k)+cos((pi/64)*(2*(k-1)+1)*(q-1-16))*cq(q);
        end
    end
    s=sk;
end