#include <openssl/pem.h>
#include <openssl/evp.h>
#include <openssl/obj_mac.h>

void handleErrors() { abort(); }

int main() {
  size_t secret_len_v;
  size_t* secret_len = &secret_len_v;

  EVP_PKEY_CTX *ctx;
  unsigned char *secret;
  EVP_PKEY *pkey = NULL, *peerkey;
  if (NULL == PEM_read_PrivateKey(fopen("alice-ec.pem", "r"), &pkey, NULL, NULL)) handleErrors();
  if (NULL == PEM_read_PUBKEY(fopen("bob-ec.pub", "r"), &peerkey, NULL, NULL)) handleErrors();
  if(NULL == (ctx = EVP_PKEY_CTX_new(pkey, NULL))) handleErrors();
  if(1 != EVP_PKEY_derive_init(ctx)) handleErrors();
  if(1 != EVP_PKEY_derive_set_peer(ctx, peerkey)) handleErrors();

  /* Determine buffer length for shared secret */
  if(1 != EVP_PKEY_derive(ctx, NULL, secret_len)) handleErrors();

  /* Create the buffer */
  if(NULL == (secret = OPENSSL_malloc(*secret_len))) handleErrors();

  /* Derive the shared secret */
  if(1 != (EVP_PKEY_derive(ctx, secret, secret_len))) handleErrors();

  fwrite(secret, secret_len_v, 1, stdout);

  EVP_PKEY_CTX_free(ctx);
  EVP_PKEY_free(peerkey);
  EVP_PKEY_free(pkey);

  return 0;
}
