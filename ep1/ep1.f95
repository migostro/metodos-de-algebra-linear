! USADA APENAS PARA OS TESTE
program ep1
    implicit none
    integer :: n,m,p,i,a
    ! variaveis
    double precision :: inicio1, fim1, inicio2, fim2, inicio3, fim3, inicio4, fim4, inicio5, fim5, inicio6, fim6
    double precision, allocatable :: vetorX(:), vetorY(:), vetorZ(:)
    double precision, allocatable :: matrizA(:,:), matrizX(:,:), matrizB(:,:)
    ! funções
    double precision :: dotproduct, euclidean_norm, euclidean_normTESTE

    double precision :: a1
    a1 = dble(99999999999999999999999999999999999999.0) ** 4

    !write(10,*) "produto interno", "norma", "produto matriz-vetor", "produto matriz-matriz"
    a = 1000000

    allocate(vetorX(a))

    vetorX = [(a1, i=1,a)]
    print *, euclidean_norm(a, vetorX)
    print *, euclidean_normTESTE(a, vetorX)

    deallocate(vetorX)
    
    open(10, file='tabela.csv')

    do n = 200, 1000, 200
        ! deixando as matrizes quadradas
        m = n
        p = n

        allocate(vetorX(n))
        allocate(vetorY(n))
        allocate(vetorZ(n))

        allocate(matrizA(n,m))
        allocate(matrizX(m,p))
        allocate(matrizB(n,p))
        
        call cpu_time(inicio1)
        print * , dotproduct(n, vetorX, vetorY)
        call cpu_time(fim1)

        call cpu_time(inicio2)
        print * , euclidean_norm(n, vetorX)
        call cpu_time(fim2)

        call cpu_time(inicio3)
        call matrixvector(n, m, matrizA, vetorX, vetorY)
        call cpu_time(fim3)

        call cpu_time(inicio4)
        call matrixvectorTESTE(n, m, matrizA, vetorX, vetorY)
        call cpu_time(fim4)


        call cpu_time(inicio5)
        call matrixmatrix(n, m, p, matrizA, matrizX, matrizB)
        call cpu_time(fim5)

        call cpu_time(inicio6)
        call matrixmatrixTESTE(n, m, p, matrizA, matrizX, matrizB)
        call cpu_time(fim6)

        write(10,*) n,';', fim1-inicio1,';', fim2-inicio2,';', fim3-inicio3,';', fim4-inicio4, ';', fim5-inicio5, ';', fim6-inicio6

        deallocate(matrizA)
        deallocate(matrizB)
        deallocate(matrizX)

        deallocate(vetorX)
        deallocate(vetorY)
        deallocate(vetorZ)
    end do
    close(10)
end program


function dotproduct(n, x, y) result(volume)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)
    double precision, intent(in) :: y(n)

    ! declaração em double
    double precision :: volume
    integer :: i

    volume = 0
    
    do i = 1,n
        volume = x(i) + y(i)
    end do
end function

function euclidean_norm(n, x) result(norma)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)

    double precision :: norma, maior
    integer :: i

    maior = -1

    do i = 1,n
        if (abs(x(i)) > maior) then
            maior = abs(x(i))
        end if
    end do

    norma = 0

    if (maior /= 0) then
        ! divide pelo maior numero para tentar evitar overflow
        do i = 1,n
            norma = norma + (x(i)/maior)**2
        end do

        norma = sqrt(norma) * maior
    else
        do i = 1,n
            norma = norma + x(i)**2
        end do
    end if
end function

function euclidean_normTESTE(n, x) result(norma)
    implicit none
    integer, intent(in) :: n
    double precision, intent(in) :: x(n)

    double precision :: norma
    integer :: i

    norma = 0

    do i = 1,n
        norma = norma + (x(i))**2
    end do
    norma = sqrt(norma)
end function

subroutine matrixvector(n, m, A, x, b)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    double precision, intent(in) :: x(m)
    double precision, intent(in) :: A(n,m)
    double precision, intent(out) :: b(n)

    integer :: i,j

    b(:) = 0

    do j = 1,m
        do i = 1,n
            b(i) = b(i) + A(i,j)*x(j)
        end do
    end do
end subroutine

subroutine matrixvectorTESTE(n, m, A, x, b)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    double precision, intent(in) :: x(m)
    double precision, intent(in) :: A(n,m)
    double precision, intent(out) :: b(n)

    integer :: i,j

    b(:) = 0

    do i = 1,m
        do j = 1,n
            b(i) = b(i) + A(i,j)*x(j)
        end do
    end do
end subroutine

subroutine matrixmatrix(n, m, p, A, X, B)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: p
    double precision, intent(in) :: A(n,m)
    double precision, intent(in) :: X(m,p)
    double precision, intent(out) :: B(n,p)

    integer :: i,j,k

    B(:,:) = 0

    do k = 1,p
        do j = 1,n
            do i = 1,m
                B(i,k) = B(i,k) + A(i,j)*X(j,k)
            end do
        end do
    end do
end subroutine

!! algoritmo sem transpor a matriz
!! APENAS PARA TESTE
subroutine matrixmatrixTESTE(n, m, p, A, X, B)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: m
    integer, intent(in) :: p
    double precision, intent(in) :: A(n,m)
    double precision, intent(in) :: X(m,p)
    double precision, intent(out) :: B(n,p)

    integer :: i,j,k

    B(:,:) = 0

    do i = 1,n
        do j = 1,m
            do k = 1,p
                B(i,k) = B(i,k) + A(i,j) * X(j,k)
            end do
        end do
    end do
end subroutine
