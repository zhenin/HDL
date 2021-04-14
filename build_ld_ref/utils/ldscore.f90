! computing LDscore using banded R info from plink
! M = number of markers
! nb = half-bandwidth above the diagonal
! plkvec = padded plink output in vector format
!
     SUBROUTINE ldscore(M,nb,plkvec,score)
     double precision plkvec((M-1)*nb), score(M)

! do this in pieces: 
! easiest: rectangular part: diag and upper diag 
     idx= 0
     score(M) = 1.
     do 10 i=1,(M-1)
       score(i) = 1.
       do 20 j=1,nb
         idx = idx+1
         score(i) = score(i) + plkvec(idx)**2
20     continue
10 continue

! lower diag: top-tail, 
     do 100 i=2,nb
       do 200 j=1,(i-1)
         idx = (i-1) +(j-1)*(nb-1)
         score(i) = score(i) + plkvec(idx)**2
200    continue
100  continue
! lower diag: main part
     do 500 i=(nb+1),M
       do 600 j=0,(nb-1)
        jx = (i-nb)+j
        idx = (i-1) + (i-nb-1+j)*(nb-1)
        score(i) = score(i) + plkvec(idx)**2
600    continue
500  continue

     RETURN
     END


    
