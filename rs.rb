class ZeroDivisionError
end

class GaloisField
	def initialize()
		setupTables()
	end
	
	def setupTables()
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

		ret = @gf_exp[@gf_log[x] + @gf_log[y]]
		#puts("gf_mul: x=0x#{x.to_s(16)}, y=0x#{y.to_s(16)}, ret=0x#{ret.to_s(16)}")
		ret
	end

	def divide(x, y)
		raise ZeroDivisionError.new if(y == 0)
		return 0 if(x == 0)

		@gf_exp[@gf_log[x] + 255 - @gf_log[y]]
	end

	def polynomialScale(p, x)
		p.map() { |i| multiply(i, x) }
	end

	def polynomialAdd(p, q)
		r = [0] * [p.size, q.size].max

		p.size.times() do |i|
			r[i + r.size - p.size] = p[i]
		end

		q.size.times() do |i|
			r[i + r.size - q.size] ^= q[i]
		end

		r
	end

	def polynomialMultiply(p, q)
		r = [0] * (p.size + q.size - 1)
		#puts(p.inspect)
		#puts(q.inspect)
		q.size.times() do |j|
			p.size.times() do |i|
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
			tomult = [1, @gf_exp[i]]
			g = polynomialMultiply(g, tomult)
		end

		g
	end

	def rsEncode(msg, nsym)
		gen = polynomialGenerator(nsym)
		msgOut = [0] * (msg.size + nsym)

		msg.size.times() do |i|
			msgOut[i] = msg[i]
		end

		msg.size.times() do |i|
			coefficient = msgOut[i]

			gen.size.times() do |j|
				msgOut[j + i] ^= multiply(gen[j], coefficient)
			end unless(coefficient == 0)
		end

		msg.size.times() do |i|
			msgOut[i] = msg[i]
		end

		#debugHexDump(msgOut, "msgOut")

		msgOut
	end

	def calcSyndromes(msg, nsym)
		synd = [0] * nsym
		
		nsym.times() do |i| 
			#puts("msg=#{msg.inspect}")
			ret = polynomialEval(msg, @gf_exp[i]) 
			#puts("ret=%02x, gf_exp[i]=%02x, i=%d" % [ret, @gf_exp[i], i])
			synd[i] = ret
		end

		return synd
	end

	def correctErrata(msg, synd, pos)
		# calculate error locator polynomial
		q = [1]
		pos.size.times() do |i|
			x = @gf_exp[msg.size - 1 - pos[i]]
			q = polynomialMultiply(q, [x, 1])
		end

		#puts("q: #{q.inspect}")

		# calculate error evaluator polynomial
		p = synd[0 .. pos.size - 1]
		#puts("selected syndromes: #{p.inspect}")
		p.reverse!()
		p = polynomialMultiply(p, q)
		#puts("p: (size=#{p.size}) #{p.inspect}")
		p = p[p.size - pos.size .. p.size - 1]
		#puts("p2: (size=#{p.size}) #{p.inspect}")
		#puts("q: #{q.inspect}")
		# eliminate even terms
		qsize = q.size
		#puts("q[#{qsize&1} .. #{qsize - 1}]")
		q = q[qsize & 1 .. qsize - 1]
		q1 = []
		q.each_index() do |i|
			q1 << q[i] if(i & 1 == 0)
		end
		q = q1

		#puts("q after every-second filtering: #{qsize & 1} #{q.inspect}")
		#
		# compute corrections
		pos.size.times() do |i|
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
				#puts("synd[%d]=%02x" % [ i - j, synd[i - j] ])
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
		if(errs * 2 > synd.size)
			# Too many errors to correct
			return nil
		end

		# Find zeros of error polynomial
		errPos = []
		nmess.times() do |i|
			if(polynomialEval(errPoly, @gf_exp[255 - i]) == 0)
				errPos << nmess - 1 - i
			end
		end

		if(errPos.size != errs)
			# can't find error locations (?)
			#puts("errPos.size=#{errPos.size}, errs=#{errs}")
			return nil
		end

		return errPos
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
		return fsynd
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

		if(erasePos.size > nsym)
			# too many erasures
			puts("too many erasures")
			return nil
		end

		synd = calcSyndromes(msgOut, nsym)
		if(synd.max() == 0)
			# no errors detected?
			#puts("no errors detected")
			return msgOut
		end

		fsynd = forneySyndromes(synd, erasePos, msgOut.size)
		errPos = findErrors(fsynd, msgOut.size)
		if(errPos == nil)
			# error location failed
			#puts("error location failed")
			return nil
		end

		correctErrata(msgOut, synd, erasePos + errPos)
		synd = calcSyndromes(msgOut, nsym)
		if(synd.max() > 0)
			# message still not right
			#puts("message still not right: synd.size = #{synd.size}")
			return nil
		end

		return msgOut
	end
end
