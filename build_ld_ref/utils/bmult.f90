! computing R%*%x using banded R info from plink
! to simplify: pad some zeroes at the end of
! -- plink output to make it rectangular 
! -- x vector
!
! M = number of markers
! nb = half-bandwidth above the diagonal
! plkvec = padded plink output in vector format
!
     SUBROUTINE bmult(M,nb,plkvec,x,out)
     double precision plkvec((M-1)*nb),x(M+nb),out(M)

! do this in pieces: 
! easiest: rectangular part: diag and upper diag 
     idx= 0
     out(M) = x(M)
     do 10 i=1,(M-1)
       out(i) = x(i)
       do 20 j=1,nb
         idx = idx+1
         out(i) = out(i) + plkvec(idx)*x(i+j)
20     continue
10 continue

! lower diag: top-tail, 
     do 100 i=2,nb
       do 200 j=1,(i-1)
         idx = (i-1) +(j-1)*(nb-1)
         out(i) = out(i) + plkvec(idx)*x(j)
200    continue
100  continue
! lower diag: main part
     do 500 i=(nb+1),M
       do 600 j=0,(nb-1)
        jx = (i-nb)+j
        idx = (i-1) + (i-nb-1+j)*(nb-1)
        out(i) = out(i) + plkvec(idx)*x(jx)
600    continue
500  continue

    RETURN
    END

