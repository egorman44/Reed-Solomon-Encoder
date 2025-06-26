package rs

import chisel3._
import chisel3.util._

class GfFieldMath(val symbWidth: Int, val primitivePoly: Int, val primitiveElement: Int = 2) {

  val elementsNum = 1 << symbWidth
  val maxExponent = elementsNum - 1

  /////////////////////////////
  // alphaToSymb maps powers of the primitive element ("alpha").
  // There are maxExponent number of non zero elements in GF that
  // could be represented as power of alpha. But use all elementsNum
  // in the table to simplify HW.
  // For example, if primitiveElement=2 and symbWidth=8:
  // 
  // tbl[0] = alpha ^ 0 = 1
  // tbl[1] = alpha ^ 1 = 2
  // ...
  // tbl[255] = 0 is invalid and should be excluded by the function
  //            that calls alphaToSymb()
  /////////////////////////////

  val alphaToSymb: Array[Int] = {
    val tbl = Array.fill(elementsNum)(0)
    tbl(0) = 1
    for (i <- 1 until maxExponent) {
      val next = tbl(i - 1) * primitiveElement
      tbl(i) = if (next >= elementsNum) next ^ primitivePoly else next
    }
    tbl
  }

  val symbToAlpha: Array[Int] = {
    val tbl = Array.fill(elementsNum)(0)
    for (i <- 0 until maxExponent) {
      tbl(alphaToSymb(i)) = i
    }
    tbl
  }

  def toSymbol(alpha: Int): Int = alphaToSymb(alpha)

  def toAlpha(symb: Int): Int = symbToAlpha(symb)

  def gfMult(a: Int, b: Int): Int = {
    if (a == 0 || b == 0) 0
    else {
      val alphaA = toAlpha(a)
      val alphaB = toAlpha(b)
      toSymbol((alphaA + alphaB) % maxExponent)
    }
  }

  def gfDiv(dividend: Int, divider: Int): Int = {
    if (dividend == 0) 0
    else if (divider == 0) throw new ArithmeticException("Division by zero in GF is undefined")
    else {
      val alphaDividend = toAlpha(dividend)
      val alphaDivider = toAlpha(divider)
      toSymbol((alphaDividend - alphaDivider + maxExponent) % maxExponent)
    }
  }

  def gfInv(symb: Int): Int = {
    if (symb == 0)
      throw new ArithmeticException("Inverse of zero is undefined in GF")
    else {
      val alpha = toAlpha(symb)
      toSymbol((maxExponent - alpha) % maxExponent)
    }
  }

  def gfPow(symb: Int, degree: Int): Int = {
    if (symb == 0) 0
    else {
      val alpha = toAlpha(symb)
      toSymbol((alpha * degree) % maxExponent)
    }
  }

  def polyMult(p: Seq[Int], q: Seq[Int]): Seq[Int] = {
    val result = Array.fill(p.length + q.length - 1)(0)
    for (i <- p.indices; j <- q.indices) {
      result(i + j) ^= gfMult(p(i), q(j))
    }
    result
  }

  def genGeneratorPoly(nsym: Int, fcr: Int = 1): Seq[Int] = {
    var gen = Seq(1)
    for (i <- 0 until nsym) {
      val factor = Seq(1, gfPow(primitiveElement, i + fcr))
      gen = polyMult(gen, factor)
    }
    gen
  }

  def computeGammaMatrix(I: Int, T: Int): Array[Array[Int]] = {
    val g = genGeneratorPoly(2 * T).reverse
    val gamma = Array.fill(2 * T, I + 1)(0)

    for (i <- 0 until 2 * T) {
      for (omega <- 1 to I) {
        val term1 = if (i - omega + 1 >= 0 && i - omega + 1 < g.length) g(i - omega + 1) else 0
        var acc = 0
        for (j <- 1 to omega) {
          val gIdx = 2 * T - j
          val gVal = if (gIdx >= 0 && gIdx < g.length) g(gIdx) else 0
          val gammaPrev = gamma(i)(omega - j)
          acc ^= gfMult(gVal, gammaPrev)
        }
        gamma(i)(omega) = term1 ^ acc
      }
    }

    gamma
  }
}
