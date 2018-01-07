
    !-------------------------------------------------------------------------------------------------------------
    !
    !> \file    GetElementName.f90
    !> \brief   Get the name of each chemical element.
    !> \author  M.H.A. Piro
    !> \date    Apr. 24, 2012
    !> \sa      CheckSystem.f90
    !
    !
    ! Revisions:
    ! ==========
    !
    !    Date          Programmer        Description of change
    !    ----          ----------        ---------------------
    !    10/18/2011    M.H.A. Piro       Original code
    !
    !
    ! Purpose:
    ! ========
    !
    !> \details The purpose of this subroutine is to return the name of each chemical element on the periodic
    !! table.
    !
    !-------------------------------------------------------------------------------------------------------------

subroutine GetElementName(cElementNamePT)

    USE ModuleThermo

    implicit none

    character(3),dimension(0:nElementsPT):: cElementNamePT


    ! Define element name (short-hand):
    cElementNamePT(0)   = 'e'
    cElementNamePT(1)   = 'H'
    cElementNamePT(2)   = 'He'
    cElementNamePT(3)   = 'Li'
    cElementNamePT(4)   = 'Be'
    cElementNamePT(5)   = 'B'
    cElementNamePT(6)   = 'C'
    cElementNamePT(7)   = 'N'
    cElementNamePT(8)   = 'O'
    cElementNamePT(9)   = 'F'
    cElementNamePT(10)  = 'Ne'
    cElementNamePT(11)  = 'Na'
    cElementNamePT(12)  = 'Mg'
    cElementNamePT(13)  = 'Al'
    cElementNamePT(14)  = 'Si'
    cElementNamePT(15)  = 'P'
    cElementNamePT(16)  = 'S'
    cElementNamePT(17)  = 'Cl'
    cElementNamePT(18)  = 'Ar'
    cElementNamePT(19)  = 'K'
    cElementNamePT(20)  = 'Ca'
    cElementNamePT(21)  = 'Sc'
    cElementNamePT(22)  = 'Ti'
    cElementNamePT(23)  = 'V'
    cElementNamePT(24)  = 'Cr'
    cElementNamePT(25)  = 'Mn'
    cElementNamePT(26)  = 'Fe'
    cElementNamePT(27)  = 'Co'
    cElementNamePT(28)  = 'Ni'
    cElementNamePT(29)  = 'Cu'
    cElementNamePT(30)  = 'Zn'
    cElementNamePT(31)  = 'Ga'
    cElementNamePT(32)  = 'Ge'
    cElementNamePT(33)  = 'As'
    cElementNamePT(34)  = 'Se'
    cElementNamePT(35)  = 'Br'
    cElementNamePT(36)  = 'Kr'
    cElementNamePT(37)  = 'Rb'
    cElementNamePT(38)  = 'Sr'
    cElementNamePT(39)  = 'Y'
    cElementNamePT(40)  = 'Zr'
    cElementNamePT(41)  = 'Nb'
    cElementNamePT(42)  = 'Mo'
    cElementNamePT(43)  = 'Tc'
    cElementNamePT(44)  = 'Ru'
    cElementNamePT(45)  = 'Rh'
    cElementNamePT(46)  = 'Pd'
    cElementNamePT(47)  = 'Ag'
    cElementNamePT(48)  = 'Cd'
    cElementNamePT(49)  = 'In'
    cElementNamePT(50)  = 'Sn'
    cElementNamePT(51)  = 'Sb'
    cElementNamePT(52)  = 'Te'
    cElementNamePT(53)  = 'I'
    cElementNamePT(54)  = 'Xe'
    cElementNamePT(55)  = 'Cs'
    cElementNamePT(56)  = 'Ba'
    cElementNamePT(57)  = 'La'
    cElementNamePT(58)  = 'Ce'
    cElementNamePT(59)  = 'Pr'
    cElementNamePT(60)  = 'Nd'
    cElementNamePT(61)  = 'Pm'
    cElementNamePT(62)  = 'Sm'
    cElementNamePT(63)  = 'Eu'
    cElementNamePT(64)  = 'Gd'
    cElementNamePT(65)  = 'Tb'
    cElementNamePT(66)  = 'Dy'
    cElementNamePT(67)  = 'Ho'
    cElementNamePT(68)  = 'Er'
    cElementNamePT(69)  = 'Tm'
    cElementNamePT(70)  = 'Yb'
    cElementNamePT(71)  = 'Lu'
    cElementNamePT(72)  = 'Hf'
    cElementNamePT(73)  = 'Ta'
    cElementNamePT(74)  = 'W'
    cElementNamePT(75)  = 'Re'
    cElementNamePT(76)  = 'Os'
    cElementNamePT(77)  = 'Ir'
    cElementNamePT(78)  = 'Pt'
    cElementNamePT(79)  = 'Au'
    cElementNamePT(80)  = 'Hg'
    cElementNamePT(81)  = 'Tl'
    cElementNamePT(82)  = 'Pb'
    cElementNamePT(83)  = 'Bi'
    cElementNamePT(84)  = 'Po'
    cElementNamePT(85)  = 'At'
    cElementNamePT(86)  = 'Rn'
    cElementNamePT(87)  = 'Fr'
    cElementNamePT(88)  = 'Ra'
    cElementNamePT(89)  = 'Ac'
    cElementNamePT(90)  = 'Th'
    cElementNamePT(91)  = 'Pa'
    cElementNamePT(92)  = 'U'
    cElementNamePT(93)  = 'Np'
    cElementNamePT(94)  = 'Pu'
    cElementNamePT(95)  = 'Am'
    cElementNamePT(96)  = 'Cm'
    cElementNamePT(97)  = 'Bk'
    cElementNamePT(98)  = 'Cf'
    cElementNamePT(99)  = 'Es'
    cElementNamePT(100) = 'Fm'
    cElementNamePT(101) = 'Md'
    cElementNamePT(102) = 'No'
    cElementNamePT(103) = 'Lr'
    cElementNamePT(104) = 'Rf'
    cElementNamePT(105) = 'Db'
    cElementNamePT(106) = 'Sg'
    cElementNamePT(107) = 'Bh'
    cElementNamePT(108) = 'Hs'
    cElementNamePT(109) = 'Mt'
    cElementNamePT(110) = 'Ds'
    cElementNamePT(111) = 'Rg'
    cElementNamePT(112) = 'Cn'
    cElementNamePT(113) = 'Uut'
    cElementNamePT(114) = 'Uuq'
    cElementNamePT(115) = 'Uup'
    cElementNamePT(116) = 'Uuh'
    cElementNamePT(117) = 'Uus'
    cElementNamePT(118) = 'Uuo'

    return

end subroutine GetElementName

