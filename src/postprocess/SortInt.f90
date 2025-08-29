subroutine SortInt(iArray, n)
    INTEGER :: i, j, k, temp, N
    integer,dimension(N):: iArray

   ! Selection Sort
    DO i = 1, N - 1
      k = i
      DO j = i + 1, N
        IF (iArray(j) < iArray(k)) THEN
          k = j
        END IF
      END DO
      IF (k /= i) THEN
        temp = iArray(i)
        iArray(i) = iArray(k)
        iArray(k) = temp
      END IF
    END DO

end subroutine SortInt