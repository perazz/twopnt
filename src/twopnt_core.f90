module twopnt
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, twopnt!"
  end subroutine say_hello
end module twopnt
