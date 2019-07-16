module m_link
  use m_shell
  implicit none

  type :: link
  class(shell), pointer :: variable => null()
  type(link), pointer :: next => null()
  end type

end module