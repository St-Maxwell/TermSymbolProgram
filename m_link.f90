module m_link
  implicit none

  type :: link
  class(*), pointer :: variable => null()
  type(link), pointer :: next => null()
  end type

end module