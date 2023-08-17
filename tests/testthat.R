library(testthat)
library(MiscMetabar)

test_check("MiscMetabar")

# expect_equal() is equal within small numerical tolerance?
# expect_identical() is exactly equal?
# expect_match() matches specified string or regular
# expect_output() prints specified output?
# expect_message() displays specified message?
# expect_warning() displays specified warning?
# expect_error() throws specified error?
# expect_is() output inherits from certain class?
# expect_false() returns FALSE?
# expect_true() returns TRUE?

##   expect_error(1 / "a", "non-numeric argument")
##   expect_warning(log(-1), "NaNs produced")