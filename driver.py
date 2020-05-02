import ReedSolomon
import fieldmath
import random
import numpy as np


def char_to_field(c, field):
    """ In the future I might add support for multiple different fields. But for now only BinaryField(ox11d) is allowed
    """
    if isinstance(c, str):
        c = ord(c[0])
    if not isinstance(c, int) or not isinstance(c, np.uint8):
        c = int(c)
    if not isinstance(field, fieldmath.BinaryField) or field.mod != 0x11d:
        raise Exception("Unsupported field")

    return c


def field_to_char(x, field):
    """ Inverse of the char_to_field function
    """
    if not isinstance(field, fieldmath.BinaryField) or field.mod != 0x11d:
        raise Exception("Unsupported field")
    x = field.multiply(field.one(), x)  # Ensures the provided value is in the field
    return chr(x)


def msg_to_field(s, field):
    return [char_to_field(c, field) for c in s]


def field_to_msg(x_arr, field):
    return ''.join([str(field_to_char(x, field)) for x in x_arr])


def split(lst, size):
    res = []
    for i in range(0, len(lst), size):
        res.append(lst[i:i + size])
    if len(res[-1]) < size:
        res[-1] = res[-1] + ' '*(size-len(res[-1]))
    return res


def sim_error(msg, field, p_err=float(1/10), p_con=0.0):
    res = [None]*len(msg)
    for m in range(len(res)):
        if m == 0 or msg[m-1] == res[m-1]:
            if random.random() < p_err:
                e = random.randint(0, field.size-1)
                res[m] = field.add(msg[m], e)
            else:
                res[m] = msg[m]
        elif msg[m-1] != res[m-1]:
            if random.random() < p_err+p_con:
                e = random.randint(0, field.size-1)
                res[m] = field.add(msg[m], e)
            else:
                res[m] = msg[m]
        else:
            res[m] = msg[m]

    return res


if __name__ == '__main__':
    f = fieldmath.BinaryField(0x11d)
    # GRS parameters
    # getting .0259
    k = 24
    n = 30

    # Random error probabilities
    # Probability of byte_error
    p_e = float(1/30)
    # Probability of consecutive_byte_error
    p_c = float(1/18)

    grs = ReedSolomon.GeneralizedReedSolomon(f=f, k=k, n=n, alpha=0x2, v_arr=1, conventional_creation=True)
    in_msg = input("Enter your message: ")
    print(f"Input message: {in_msg}")

    msgs = split(in_msg, k)
    in_msg = "".join(msgs)

    encoded_msg = []
    for j in range(len(msgs)):
        encoded_msg += grs.encode(msg_to_field(msgs[j], f))

    print(f"Encoded message   : {encoded_msg}")

    # Simulate Error
    errors = 0
    total = 1000
    for i in range(total):
        encoded_msg_with_error = sim_error(encoded_msg, f, p_err=p_e, p_con=p_c)
        error = [f.subtract(encoded_msg_with_error[i], encoded_msg[i]) for i in range(len(encoded_msg))]
        encoded_msg_with_error = split(encoded_msg_with_error, n)

        decoded_msgs = []
        for j in range(len(encoded_msg_with_error)):
            decoded_msgs += grs.decode(encoded_msg_with_error[j])

        decoded_msg_str = field_to_msg(decoded_msgs, f)
        if ''.join(decoded_msg_str) == in_msg:
            print(f"{i} Successful decoding: True")
        else:
            print(f"{i} Successful decoding: False")
            print(f"Error:{error}")

            errors += 1
    print(f"Error rate: {float(errors/total)}")
