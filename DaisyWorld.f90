program daisyworld
    implicit none
    real, parameter :: aw=0.75, ab=0.25, ag=0.5, gamma=0.3, p=1.0, q=20.0, sigma=5.670374419e-8
    real :: L=0.50, S=917.0
    real :: albedo
    real :: alpha_w=0.01, alpha_b=0.01, alpha_g
    real :: alpha_wc, alpha_bc, alpha_wold, alpha_bold
    real :: beta_w, beta_b
    real :: temp_w, temp_b, temp_e
    real :: delalpha_w, delalpha_b, deltime=1.0
    integer :: i, imax=1000
    real, parameter :: tollerance=0.000001
    
    open(unit=100, file="daisyworld.txt", status='replace', action='write')
    write(100,'(5A10)') "sunLum","tPlanet","aWhite","aBlack","aGround"
    do while(L <= 1.65)
        if(alpha_w < 0.01) alpha_w=0.01
        if(alpha_b < 0.01) alpha_b=0.01
        alpha_g = p - alpha_w - alpha_b
        alpha_wc = tollerance
        alpha_bc = tollerance
        i = 0
        do while ((i <= imax) .and. (alpha_wc >= tollerance) .and. (alpha_bc >= tollerance))
           call CAL_ALBEDO(alpha_w, alpha_b, alpha_g, aw, ab, ag, albedo)
           call CAL_TEMP_E(S, L, albedo, sigma, temp_e) 
           call CAL_TEMP(albedo, aw, temp_e, q, temp_w)
           call CAL_TEMP(albedo, ab, temp_e, q, temp_b)
           call CAL_BETA(temp_w, beta_w)
           call CAL_BETA(temp_b, beta_b)
           call CAL_DelALPHA(alpha_w, alpha_g, beta_w, gamma, delalpha_w)
           call CAL_DelALPHA(alpha_b, alpha_g, beta_b, gamma, delalpha_b)
           alpha_wold = alpha_w
           alpha_bold = alpha_b
           alpha_w = alpha_w + delalpha_w*deltime
           alpha_b = alpha_b + delalpha_b*deltime
           alpha_g = p - alpha_w - alpha_b
           alpha_wc = abs(alpha_w - alpha_wold)
           alpha_bc = abs(alpha_b - alpha_bold)
           i = i + 1
        end do
        write(100,'(5F10.5)') L,temp_e-273,alpha_w*100,alpha_b*100,alpha_g*100
        L = L + 0.001
    end do
end program daisyworld

subroutine CAL_ALBEDO(alphaw, alphab, alphag, aw, ab, ag, albedo)
    real, intent(in) :: alphaw, alphab, alphag, aw, ab, ag
    real, intent(out) :: albedo
    albedo = alphaw*aw + alphab*ab + alphag*ag
end subroutine CAL_ALBEDO

subroutine CAL_DelALPHA(alpha, alphag, beta, gamma, delalpha)
    real, intent(in) :: alpha, alphag, beta, gamma
    real, intent(out) :: delalpha
    delalpha = alpha*(alphag*beta - gamma)
end subroutine CAL_DelALPHA

subroutine CAL_BETA(temp, beta)
    real, intent(in) :: temp
    real, intent(out) :: beta
    if (temp > 278.0 .and. temp < 313.0) then
        beta = 1.0 - 0.003265*(295.5-temp)**2.0
    else
        beta = 0.0
    end if 
end subroutine CAL_BETA

subroutine CAL_TEMP_E(S, L, albedo_a, sigma, temp_e)
    real, intent(in) :: S, L, albedo_a, sigma
    real, intent(out) :: temp_e
    temp_e = ((S*L*(1-albedo_a))/sigma)**(1.0/4.0)
end subroutine CAL_TEMP_E

subroutine CAL_TEMP(albedo_a, albedo, temp_e, q, temp)
    real, intent(in) :: albedo_a, albedo, temp_e, q
    real, intent(out) :: temp
    temp = q*(albedo_a-albedo) + temp_e
end subroutine CAL_TEMP
