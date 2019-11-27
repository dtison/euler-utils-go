package utils

import (
	"container/list"
	"fmt"
	"github.com/yourbasic/bit"
	"math"
	"math/big"
	"reflect"
	"strconv"
	"time"
)

//const Pi = 3.14
var int64Type = reflect.TypeOf(int64(0))

func makeInt64(val interface{}) int64 {
	return reflect.ValueOf(val).Convert(int64Type).Int()
}

// GetDivisors creates a list of divisors for number
func GetDivisors(number int) []int {
	var result, needReverse []int
	if number == 1 {
		result = append(result, 1)
	}
	for i := 1; i <= int(math.Sqrt(float64(number))); i++ {
		if number%i == 0 {
			result = append(result, i)
			quotient := number / i
			if quotient != i {
				needReverse = append(needReverse, quotient)
			}
		}
	}
	for i := range needReverse {
		result = append(result, needReverse[len(needReverse)-i-1])
	}
	return result
}

// GetProperDivisors creates a list/slice of properdivisors for number
func GetProperDivisors(number int) []int {
	divisors := GetDivisors(number)
	//	fmt.Println("Proper Divisors of ", number, divisors[0:len(divisors)-1])
	return divisors[0 : len(divisors)-1]
}

// ArraySum returns sum of all numbers in array
func ArraySum(array []int) int {
	sum := 0
	for _, v := range array {
		sum += v
	}
	return sum
}

/* func GetPrimeFactors(number interface{}) (*list.List, bool) {
	n := makeInt64(number)
	results := list.New()
	for i := int64(2); i <= n/i; i++ {
		for n%i == 0 {
			results.PushBack(i)
			n /= i
		}
	}
	if n > 1 {
		results.PushBack(n)
	}
	return results, results.Len() > 0
} */

// GetPrimeFactors creates a list of prime factors for number
func GetPrimeFactors(number interface{}, lastValue ...*int64) (*list.List, bool) {

	n := makeInt64(number)
	start := int64(2)
	if len(lastValue) > 0 {
		start = *lastValue[0]
	}

	results := list.New()
	value := start
	for ; value <= n/value; value++ {
		for n%value == 0 {
			results.PushBack(value)
			n /= value
		}
	}
	if n > 1 {
		results.PushBack(n)
	}
	if len(lastValue) > 0 {
		*lastValue[0] = value
	}
	return results, results.Len() > 0
}

// DisplayElapsedTime displays time since start
func DisplayElapsedTime(start time.Time) {
	fmt.Printf("Calculated in %s\n", time.Since(start))
}

// ReverseString reverses s
// based on https://github.com/golang/example/blob/master/stringutil/reverse.go
func ReverseString(s string) string {

	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < len(runes)/2; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

// FactorialFloat64 float64 factorial
func FactorialFloat64(value float64) (result float64) {
	for result = 1.0; value >= 1; value -= 1.0 {
		result = result * value
	}
	return result
}

// FactorialUint64 uint64 factorial
func FactorialUint64(value uint64) (result uint64) {
	for result = 1; value >= 1; value-- {
		result = result * value
	}
	return result
}

// Factorial int64 factorial
func Factorial(number interface{}) (result int64) {
	value := makeInt64(number)
	for result = 1; value >= 1; value-- {
		result = result * value
	}
	return result
}

// Primes means not divisible of course
type Primes interface {
	IsPrime(value int) bool
	NextPrime() int
	PrimeFactors(number interface{}) (int, *list.List)
}

// PrimeGenerator is for primes
type PrimeGenerator struct {
	limit int
	sieve *bit.Set
}

// Factor for prime factorization
type Factor struct {
	Number, Exponent int
}

// NewPrimeGenerator makes a new one
func NewPrimeGenerator(limit int) *PrimeGenerator {
	generator := new(PrimeGenerator)
	generator.limit = limit
	sieve := bit.New().AddRange(2, limit)
	for p := 2; p <= int(math.Sqrt(float64(limit))); p = sieve.Next(p) {
		for k := p * p; k < limit; k += p {
			sieve.Delete(k)
		}
	}
	generator.sieve = sieve
	return generator
}

// IsPrime returns if value is prime
func (generator *PrimeGenerator) IsPrime(value int) bool {
	return generator.sieve.Contains(value)
}

// NextPrime returns next prime
func (generator *PrimeGenerator) NextPrime(next int) int {
	return generator.sieve.Next(next)
}

// PrevPrime returns previous prime
func (generator *PrimeGenerator) PrevPrime(prev int) int {
	return generator.sieve.Prev(prev)
}

// PrimeFactors finds prime factors  TODO: Error if not enough primes
// Returns number of unique factors, and List of Factors with exponents
func (generator *PrimeGenerator) PrimeFactors(value interface{}) (int, *list.List) {

	number := makeInt64(value)
	numFactors, remainder, results := 0, number, list.New()
	if generator.sieve.Contains(int(number)) {
		return 0, results
	}

	for prime := 2; prime != -1 && remainder != 1; prime = generator.sieve.Next(prime) {

		checkFactors := func(value int64, remainder *int64) (bool, int) {
			hasFactor, exponent := false, 0
			for *remainder%value == 0 {
				*remainder /= value
				hasFactor, exponent = true, exponent+1
			}
			if hasFactor {
				numFactors++
			}
			return hasFactor, exponent
		}

		if isPrimeFactor, exponent := checkFactors(int64(prime), &remainder); isPrimeFactor {
			results.PushBack(Factor{Number: prime, Exponent: exponent})
		}
	}

	return numFactors, results
}

// IsPrime is a fallback standard prime checker
func IsPrime(n uint64) bool {
	if n <= 1 {
		return false
	}
	if n <= 3 {
		return true
	}
	if n%2 == 0 || n%3 == 0 {
		return false
	}
	for i := uint64(5); i*i <= n; i = i + 6 {
		if n%i == 0 || n%(i+2) == 0 {
			return false
		}
	}
	return true
}

// BigPow raises n to exp  TODO (Probably deprecated because of big.Int.Exp())
func BigPow(n int64, exp int64) *big.Int {
	result := big.NewInt(n)
	for i := int64(0); i < exp-1; i++ {
		result = result.Mul(result, big.NewInt(n))
	}
	return result
}

// TODO: Public Remove element from slice
func remove(slice []int, s int) []int {
	return append(slice[:s], slice[s+1:]...)
}

// GetDigits returns digits
func GetDigits(source interface{}) (result []byte) {
	number := makeInt64(source)
	for number > 0 {
		result = append(result, byte(number%10))
		number /= 10
	}
	// Reverse the array
	for left, right := 0, len(result)-1; left < right; left, right = left+1, right-1 {
		result[left], result[right] = result[right], result[left]
	}
	return
}

// GetDigitsASCII same as GetDigits but converts them all to ascii
func GetDigitsASCII(source interface{}) (result []byte) {
	digits := GetDigits(source)
	ascii := make([]byte, len(digits))
	copy(ascii, digits)
	for i := range ascii {
		ascii[i] += 48
	}
	return ascii
}

// ArePandigitalTogether true if all parameters together include exactly 1..9
func ArePandigitalTogether(numbers ...int) bool {

	digits := new(bit.Set)

	for _, v := range numbers {
		for v > 0 {
			digit := int(v % 10)
			if digits.Contains(digit) {
				return false
			}
			digits.Add(digit)
			v /= 10
		}
	}

	if digits.Size() != 9 {
		return false
	}
	return digits.Equal(new(bit.Set).AddRange(1, 10))
}

// IsPandigital true if number has exactly 1 digit from panMin..panMax
func IsPandigital(number interface{}, panMin, panMax int) bool {

	a := makeInt64(number)
	digits := new(bit.Set)

	for a > 0 {
		digit := int(a % 10)
		if digits.Contains(digit) {
			return false
		}
		digits.Add(digit)
		a /= 10
	}

	if digits.Size() != panMax-panMin+1 {
		return false
	}
	return digits.Equal(new(bit.Set).AddRange(panMin, panMax+1))
}

// TODO:  Try using the function above here

// ArePandigitalCombined true if all three 1..9
func ArePandigitalCombined(a, b uint) bool {

	digits := new(bit.Set)
	c := a * b

	for a > 0 {
		digit := int(a % 10)
		if digits.Contains(digit) {
			return false
		}
		digits.Add(digit)
		a /= 10
	}

	for b > 0 {
		digit := int(b % 10)
		if digits.Contains(digit) {
			return false
		} else {
			digits.Add(digit)
		}
		b /= 10
	}
	for c > 0 {
		digit := int(c % 10)
		if digits.Contains(digit) {
			return false
		} else {
			digits.Add(digit)
		}
		c /= 10
	}
	if digits.Size() != 9 {
		return false
	}
	return digits.Equal(new(bit.Set).AddRange(1, 10))
}

// GeneratePrimes makes primes
func GeneratePrimes(limit int) {
	n := limit
	sieve := bit.New().AddRange(2, n)
	sqrtN := int(math.Sqrt(float64(n)))
	for p := 2; p <= sqrtN; p = sieve.Next(p) {
		for k := p * p; k < n; k += p {
			sieve.Delete(k)
		}
	}
	fmt.Println(sieve)
}

// Gcd finds greatest common divisor
func Gcd(a, b uint64) uint64 {
	x, y := a, b
	if b > a {
		y, x = a, b
	}
	for x%y != 0 {
		x, y = y, x%x
	}
	return y
}

// SwapBytesAt swaps at the indices
func SwapBytesAt(s []byte, i, j int) {
	s[i], s[j] = s[j], s[i]
}

// AppendNumbers appends 2 numbers concat
func AppendNumbers(first uint64, second uint64) uint64 {
	temp := second
	for temp > 0 {
		first *= 10
		temp /= 10
	}
	return first + second
}

// NumberDigits calculates number digits
func NumberDigits(n uint64) (result int) {
	for ; n > 0; result++ {
		n /= 10
	}
	return
}

// LeftShifts creates array of left shifts
func LeftShifts(n uint64) []uint64 {
	result := make([]uint64, 1)
	result[0] = n

	digits := NumberDigits(n)
	powTen := math.Pow(10, float64(digits-1))

	for i := 0; i < digits-1; i++ {
		firstDigit := uint64(math.Floor(float64(n) / powTen))
		left := ((n * 10) + firstDigit) -
			(firstDigit * uint64(powTen) * 10)

		result = append(result, left)
		n = left
	}
	return result
}

// IsPalindromicForBase returns true if number is palondromic in base
func IsPalindromicForBase(number uint64, base int) bool {
	reversed := uint64(0)
	for k := number; k > 0; k /= uint64(base) {
		reversed = uint64(base)*reversed + k%uint64(base)
	}
	return reversed == number
}

// GetTruncations returns array of all truncations of n, including n  Ex  3797, 797, 97, 379, 37, 7, 3.
func GetTruncations(n int) []int {

	results := make([]int, NumberDigits(uint64(n))*2-1)
	leftTruncation, rightTruncation, multiplier, index := 0, n, 1, 0

	for rightTruncation > 0 {
		if leftTruncation > 0 {
			results[index], index = leftTruncation, index+1
		}
		results[index], index = rightTruncation, index+1
		leftTruncation += multiplier * (rightTruncation % 10)
		rightTruncation /= 10
		multiplier *= 10
	}
	return results
}

// Triplet pythagorean triplet
type Triplet struct {
	A, B, C uint64
}

// PythagoreanTripletsForPerimeter returns List of triplets
func PythagoreanTripletsForPerimeter(p uint64) (*list.List, bool) {
	results := list.New()
	for a := uint64(2); a < p/3; a++ {
		term1, term2 := p*(p-2*a), 2*(p-a)
		if term1%term2 == 0 {
			triplet := new(Triplet)
			b := term1 / term2
			triplet.A, triplet.B, triplet.C = a, b, p-b-a
			results.PushBack(triplet)
		}
	}
	return results, results.Len() > 0
}

// IsTriangularNumber bool optimized from comment in https://www.mathblog.dk/project-euler-42-triangle-words/
func IsTriangularNumber(v interface{}) bool {
	n := makeInt64(v) << 1
	sqrt := int64(math.Sqrt(float64(n)))
	return n == sqrt*(sqrt+1)
}

// IsPentagonal true if is pentaganal
func IsPentagonal(number interface{}) bool {
	val := math.Sqrt(float64(makeInt64(number)*24 + 1))
	return math.Floor(val) == val && int(val)%6 == 5
}

// NthPermutationOfDigits makes //TODO: pass digits array in
func NthPermutationOfDigits(number interface{}) (result string) {
	digits := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
	n, amountRemaining := len(digits), makeInt64(number)-1

	for i := 1; i < n && amountRemaining != 0; i++ {
		factorial := Factorial(n - i)
		j := amountRemaining / factorial
		amountRemaining = amountRemaining % factorial
		result += strconv.FormatInt(int64(digits[j]), 10)
		digits = remove(digits, int(j))
	}

	for i := range digits {
		result += strconv.FormatInt(int64(digits[i]), 10)
	}
	return result
}

// Permutations calls f with each permutation of a.
func Permutations(a []rune, f func([]rune)) {
	permRecursive(a, f, 0)
}

// Permute the values at index i to len(a)-1.
func permRecursive(a []rune, f func([]rune), i int) {
	if i > len(a) {
		f(a)
		return
	}
	permRecursive(a, f, i+1)
	for j := i + 1; j < len(a); j++ {
		a[i], a[j] = a[j], a[i]
		permRecursive(a, f, i+1)
		a[i], a[j] = a[j], a[i]
	}
}

// private bool IsPalindrome(int number, int b){
// 	int reversed = 0;
// 	int k = number;

// 	while (k > 0) {
// 		reversed = b * reversed + k % b;
// 		k /= b;
// 	}

// 	bool result = (number == reversed);
// 	if (result) {
// 		Console.WriteLine("IsPalindrome {0} {1}", Convert.ToString(number, 2), reversed);
// 	}
// 	return result;
// }

/*  WIP  Permutations

func getPermutations(s []byte, l, r int) [][]byte {
	//var result [][]byte
	result := make([][]byte, 0)
	if l == r {
		fmt.Println("storing temp ans", s)
		result = append(result, s)
		return result
	} else {
		for i := l; i <= r; i++ {
			utils.SwapBytesAt(s, l, i)

			result = append(result, getPermutations(s, l+1, r)...)
			utils.SwapBytesAt(s, l, i)

		}
	}
	return result
}

*/

/* func bigFactorial(number big.Int) big.Int {
	total := new(big.Int)
	one := new(big.Int)
	one.SetString("1", 10)
	for value := number; value.Cmp(one) > 1; value.Sub(&value, one) {
		total = total.Mul(&value, total)

	}

} */

/*
  Based on this
https://jsfiddle.net/0tryqv58/

function getDivisors(n)
{
	result = [];
    // Note that this loop runs till square root
    for (let i=1; i<=Math.sqrt(n); i++)
    {
        if (n%i == 0)
        {
            // If divisors are equal, print only one
            if (n/i == i) {
            	result.push(i);
            }  else {
				result.push(i);
              result.push(n / i)
            }
        }
    }
    return result
}
*/

/*

// IsPandigital true if 1..9
func IsPandigital(n uint) bool {
	if n > 987654321 {
		return false
	}
	digits := new(bit.Set)

	for n > 0 {
		digit := int(n % 10)
		if digits.Contains(digit) {
			return false
		} else {
			digits.Add(digit)
		}
		n /= 10
	}
	return digits.Equal(new(bit.Set).AddRange(1, 10))
}

*/

/*
// Slow - tries to return in normal order
func getDigits(source uint64) (result []uint8) {
	number := source
	for number > 0 {
		// Source https://github.com/golang/go/wiki/SliceTricks
		result = append([]uint8{uint8(number % 10)}, result...)
		number /= 10
	}
	return
}

*/
/*  OLD version

func PythagoreanTriplets(limit uint64) [][]uint64 {

	result := make([][]uint64, 0)

	// triplet: a^2 + b^2 = c^2
	var a, b, c uint64
	m := uint64(2)
	c = limit
	// loop from 2 to max_limitit
	// int m = 2;

	// Limiting c would limit
	// all a, b and c
	for c <= limit {

	found:
		// now loop on j from 1 to i-1
		for n := uint64(1); n < m; n++ {

			// Evaluate and print triplets using
			// the relation between a, b and c
			a = m*m - n*n
			b = 2 * m * n
			c = m*m + n*n

			if c > limit {
				break found
			}

			fmt.Printf("%d %d %d\n", a, b, c)
		}
		m++
	}
	return result
}


*/

/*  Originally From geeks for geeks, slower

// IsTriangularNumber bool
func IsTriangularNumber1(v interface{}) bool {
	n := makeInt64(v)
	if n < 0 {
		return false
	}
	// Considering the equation n*(n+1)/2 = num
	// The equation is  : a(n^2) + bn + c = 0";
	c := (-2 * n)
	b, a := int64(1), int64(1)
	d := (b * b) - (4 * a * c)

	if d < 0 {
		return false
	}
	// Find roots of equation
	sqrt := math.Sqrt(float64(d))
	denominator := float64((2 * a))
	root1 := (float64(-b) + sqrt) / denominator
	root2 := (float64(-b) - sqrt) / denominator

	// check if roots are natural
	if root1 > 0 && math.Floor(root1) == root1 {
		return true
	}
	if root2 > 0 && math.Floor(root2) == root2 {
		return true
	}
	return false
}
*/
/*

// FactorialInt int factorial - DEPRECATED
func FactorialInt(value int) (result int) {
	for result = 1; value >= 1; value-- {
		result = result * value
	}
	return result
}

*/
