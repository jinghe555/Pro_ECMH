#先把给定的消息的hash值求出来，256bits，把它转成整数类型，作为椭圆曲线上的点的x值
#求出它对应的在椭圆曲线上的y值就得到了这个点了，消息就被映射到椭圆曲线上了
from random import choice
from gmssl import sm3,func
import cmath
#椭圆曲线方程:y^2=x^3+a*x+b
a=int('FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFC',16)
b=int('28E9FA9E9D9F5E344D5A9E4BCF6509A7F39789F515AB8F92DDBCBD414D940E93',16)
p=int('FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF00000000FFFFFFFFFFFFFFFF',16)
n=int('FFFFFFFEFFFFFFFFFFFFFFFFFFFFFFFF7203DF6B21C6052B53BBF40939D54123',16)
message='adcb652665165'
#十六进制转换为byte数组
def hex_byte(msg):
    ml=len(msg)
    if ml%2!=0:
        msg='0'+msg
    ml=int(len(msg)/2)
    msg_byte=[]
    for i in range(ml):
        msg_byte.append(int(msg[i*2:i*2+2],16))
    return msg_byte
x=int(sm3.sm3_hash(hex_byte(message)),16)
y=hex(int(pow(((pow(x,3)%p+(a*x)+b)),0.5)))
G=hex(x)[2:]+y[2:]
#随机产生字符串,OK
def random_string(strlen):
    letterlist=['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']
    strl=""
    for i in range(strlen):
        temp=choice(letterlist)
        strl+=temp
    return strl
#椭圆曲线点的计算

#倍点,OK
def doublepoint(point,length):
    l=len(point)
    len_2=2*length
    if l<len_2:
        return None
    else:
        x1=int(point[0:length],16)
        y1=int(point[length:len_2],16)
        if l==len_2:
            z1=1
        else:
            z1=int(point[len_2:],16)
        T6=(z1*z1)%p
        T2=(y1*y1)%p
        T3=(x1+T6)%p
        T4=(x1-T6)%p
        T1=(T3*T4)%p
        T3=(y1*z1)%p
        T4=(T2*8)%p
        T5=(x1*T4)%p
        T1=(T1*3)%p
        T6=(T6*T6)%p
        T6=(((a+3)%p)*T6)%p
        T1=(T1+T6)%p
        z3=(T3+T3)%p
        T3=(T1*T1)%p
        T2=(T2*T4)%p
        x3=(T3-T5)%p
        if (T5%2)==1:
            T4=(T5+((T5+p)>>1)-T3)%p
        else:
            T4=(T5+(T5>>1)-T3)%p
        T1=(T1*T4)%p
        y3=(T1-T2)%p
        form='%%0%dx'%length
        form=form*3
        return form%(x3,y3,z3)

#点加,OK
def addpoint(p1,p2,length):
    len_2=2*length
    l1=len(p1)
    l2=len(p2)
    if (l1<len_2) or (l2<len_2):
        return None
    else:
        x1=int(p1[0:length],16)
        y1=int(p1[length:len_2],16)
        if (l1==len_2):
            z1=1
        else:
            z1=int(p1[len_2:],16)
        x2=int(p2[0:length],16)
        y2=int(p2[length:len_2],16)
        T1=(z1*z1)%p
        T2=(y2*z1)%p
        T3=(x2*T1)%p
        T1=(T1*T2)%p
        T2=(T3-x1)%p
        T3=(T3+x1)%p
        T4=(T2*T2)%p
        T1=(T1-y1)%p
        z3=(z1*T2)%p
        T2=(T2*T4)%p
        T3=(T3*T4)%p
        T5=(T1*T1)%p
        T4=(x1*T4)%p
        x3=(T5-T3)%p
        T2=(y1*T2)%p
        T3=(T4-x3)%p
        T1=(T1*T3)%p
        y3=(T1-T2)%p
        form='%%0%dx'%length
        form=form*3
        return form%(x3,y3,z3)

#Jacobian加重射影坐标转换成仿射坐标,OK
def convertJacb2Nor(point,length):
    len_2=length*2
    x=int(point[0:length],16)
    y=int(point[length:len_2],16)
    z=int(point[len_2:],16)
    z_inv=pow(z,p-2,p)
    z_invSquar=(z_inv*z_inv)%p
    z_invQube=(z_invSquar*z_inv)%p
    x_new=(x*z_invSquar)%p
    y_new=(y*z_invQube)%p
    z_new=(z*z_inv)%p
    if z_new==1:
        form='%%0%dx'%length
        form=form*2
        return form%(x_new,y_new)
    else:
        print("point at infinity!!!")
        return None
    '''form='%%0%dx'%length
    form=form*2
    return form%(x_new,y_new)'''
    
#kp,点乘,OK
def kp(k,point,length):
    point='%s%s'%(point,'1')
    #print("point",point)
    mask_str='8'
    for i in range(length-1):
        mask_str+='0'
    #print("mask_str",mask_str)
    mask=int(mask_str,16)
    temp=point
    flag=False
    for n in range(length*4):
        if (flag):
            temp=doublepoint(temp,length)
        if (k&mask)!=0:
            if(flag):
                temp=addpoint(temp,point,length)
            else:
                flag=True
                temp=point
        k=k<<1
    return convertJacb2Nor(temp,length)
#签名算法
#预处理1,计算z值
def pre1(PA):
    data='001031323334353637383132333435363738'
    data+=str(a)
    data+=str(b)
    data+=G
    data+=PA
    #print("data",data)
    data_byte=hex_byte(data)
    return sm3.sm3_hash(data_byte)
#预处理2,得到杂凑值H,z：z值  m:消息
def pre2(z,m):
    data=z+m
    data_byte=hex_byte(data)
    return sm3.sm3_hash(data_byte)
#生成签名，E:消息的hash值，dA:签名者的私钥，K：随机数
def sign(E,dA,K,length,hexstr=0):
    if hexstr:
        e=int(E,16)
    else:
        E=E.encode('utf-8')
        E=E.hex()
        e=int(E,16)
    d=int(dA,16)
    k=int(K,16)
    p1=kp(k,G,length)
    x=int(p1[0:length],16)
    R=((e+x)%n)
    if R==0 or R+k==n:
        return None
    d_1=pow(d+1,n-2,n)
    S=(d_1*(k+R)-R)%n
    if S==0:
        return None
    else:
        return '%064x%064x'%(R,S)
if __name__=='__main__':
    d=random_string(64)
    k=random_string(64)
    Pa=kp(int(d,16),G,64)
    z1=pre1(Pa)
    hash_M=pre2(z1,message)
    sig=sign(hash_M,d,k,64,1)
    print("sig:",sig)
