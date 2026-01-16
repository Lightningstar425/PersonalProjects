
int bit1 = 8;
int bit2 = 11;
int co = 9;

void setup() {
 pinMode(bit1, OUTPUT);
 pinMode(bit2, OUTPUT);
 pinMode(co, OUTPUT);
}

void loop() {
  digitalWrite(bit1, HIGH);
  digitalWrite(bit2, LOW);
  digitalWrite(co, HIGH);

}
