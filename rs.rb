class ZeroDivisionError
end

class GaloisField
	def initialize()
		@gf_exp = [1] * 512
		@gf_log = [0] * 256
		x = 1

		1.upto(255) do |i|
			x <<= 1

			x = x ^ 0x11D if(x & 0x100 != 0) # in GF(2**8) subtraction is done by XOR

			@gf_exp[i] = x
			@gf_log[x] = i
		end

		255.upto(512) do |i|
			@gf_exp[i] = @gf_exp[i - 255]
		end
	end

	def multiply(x, y)
		return 0 if(x == 0 or y == 0)

		@gf_exp[@gf_log[x] + @gf_log[y]]
	end

	def divide(x, y)
		raise ZeroDivisionError.new if(y == 0)
		return 0 if(x == 0)

		@gf_exp[@gf_log[x] + 255 - @gf_log[y]]
	end

	def polynomialScale(p, x)
		p.map() do |i| 
			multiply(i, x) 
		end
	end

	def polynomialAdd(p, q)
		r = [0] * [p.size, q.size].max

		p.each_index() do |i|
			r[i + r.size - p.size] = p[i]
		end

		q.each_index() do |i|
			r[i + r.size - q.size] ^= q[i]
		end

		r
	end

	def polynomialMultiply(p, q)
		r = [0] * (p.size + q.size - 1)

		q.each_index() do |j|
			p.each_index() do |i|
				r[i + j] ^= multiply(p[i], q[j])
			end
		end

		r
	end

	# Horner scheme is used to evaluate polynomial
	def polynomialEval(p, x)
		y = p[0]

		1.upto(p.size - 1) do |i|
			y = multiply(y, x) ^ p[i]
		end

		y
	end

	def polynomialGenerator(nsym)
		g = [1]
		
		nsym.times() do |i|
			g = polynomialMultiply(g, [1, @gf_exp[i]])
		end

		g
	end

	def rsEncode(msg, nsym)
		msgOut = [0] * (msg.size + nsym)

		msg.each_index() do |i|
			msgOut[i] = msg[i]
		end

		gen = polynomialGenerator(nsym)
		msg.each_index() do |i|
			coefficient = msgOut[i]

			gen.each_index() do |j|
				msgOut[j + i] ^= multiply(gen[j], coefficient)
			end unless(coefficient == 0)
		end

		msg.each_index() do |i|
			msgOut[i] = msg[i]
		end

		msgOut
	end

	def calcSyndromes(msg, nsym)
		synd = [0] * nsym
		
		synd.each_index() do |i|
			synd[i] = polynomialEval(msg, @gf_exp[i])
		end

		synd
	end

	def correctErrata(msg, synd, pos)
		# calculate error locator polynomial
		q = [1]
		pos.each_index() do |i|
			x = @gf_exp[msg.size - 1 - pos[i]]
			q = polynomialMultiply(q, [x, 1])
		end

		# calculate error evaluator polynomial
		p = synd[0 .. pos.size - 1]
		p.reverse!()
		p = polynomialMultiply(p, q)
		p = p[p.size - pos.size .. p.size - 1]

		# eliminate even terms
		qsize = q.size

		q = q[qsize & 1 .. qsize - 1]
		q1 = []
		q.each_index() do |i|
			q1 << q[i] if(i & 1 == 0)
		end
		q = q1

		# compute corrections
		pos.each_index() do |i|
			x = @gf_exp[pos[i] + 256 - msg.size]
			y = polynomialEval(p, x)
			z = polynomialEval(q, multiply(x, x))
			msg[pos[i]] ^= divide(y, multiply(x, z))
		end
	end

	def findErrors(synd, nmess)
		# find error locator polynomial, Berlekamp-Massey algo
		errPoly = [1]
		oldPoly = [1]

		synd.each_index() do |i|
			oldPoly << 0
			delta = synd[i]

			1.upto(errPoly.size - 1) do |j|
				delta ^= multiply(errPoly[errPoly.size - 1 - j], synd[i - j])
			end

			if(delta != 0)
				if(oldPoly.size > errPoly.size)
					newPoly = polynomialScale(oldPoly, delta)
					oldPoly = polynomialScale(errPoly, divide(1, delta))
					errPoly = newPoly
				end

				errPoly = polynomialAdd(errPoly, polynomialScale(oldPoly, delta))
			end
		end

		errs = errPoly.size - 1

		# Too many errors to correct
		return nil if(errs * 2 > synd.size)

		# Find zeros of error polynomial
		errPos = []
		nmess.times() do |i|
			errPos << nmess - 1 - i if(polynomialEval(errPoly, @gf_exp[255 - i]) == 0)
		end

		# can't find error locations (?)
		return nil if(errPos.size != errs)

		errPos
	end

	def forneySyndromes(synd, pos, nmess)
		fsynd = synd.dup()

		pos.each_index() do |i|
			x = @gf_exp[nmess - 1 - pos[i]]

			(fsynd.size - 1).times() do |j|
				fsynd[j] = multiply(fsynd[j], x) ^ fsynd[j + 1]
			end

			fdynd.pop()
		end

		fsynd
	end

	def correctMessage(msgIn, nsym)
		msgOut = msgIn.dup()
		erasePos = []

		msgOut.each_index() do |i|
			if(msgOut[i] < 0)
				msgOut[i] = 0
				erasePos << i
			end
		end

		# too many erasures
		return nil if(erasePos.size > nsym)

		synd = calcSyndromes(msgOut, nsym)
		return msgOut if(synd.max() == 0) # no errors detected?

		fsynd = forneySyndromes(synd, erasePos, msgOut.size)
		errPos = findErrors(fsynd, msgOut.size)
		return nil if(errPos == nil) # error location failed

		correctErrata(msgOut, synd, erasePos + errPos)
		synd = calcSyndromes(msgOut, nsym)
		return nil if(synd.max() > 0) # message still not right

		msgOut
	end
end
